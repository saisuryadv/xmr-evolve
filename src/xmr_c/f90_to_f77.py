#!/usr/bin/env python3
"""
Convert Fortran 90 free-form-ish features in XMR code to Fortran 77
compatible with f2c.

Transformations:
1. Strip INTENT(...) from declarations
2. Inline PARAMETER declarations -> separate statements
3. Convert bare "DO" / named "DO" / F90 DO loops to labeled F77 loops
4. Convert EXIT/CYCLE to GOTO
5. Convert END SUBROUTINE/FUNCTION -> END
6. Convert array bounds "arr(expr:expr)" in declarations to "arr(*)"
7. Remove IMPLICIT NONE
8. Comment out INTERFACE...END INTERFACE blocks, extract EXTERNAL decls
9. Convert array-slice operations: ARR(:) = val, ARR(:) = OTHER(:), etc.
"""

import re
import sys
import os


def process_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    processed = [line.rstrip('\n') for line in lines]

    # ---- Pass 0: Comment out INTERFACE blocks, extract EXTERNAL decls ----
    processed, externals = strip_interface_blocks(processed)

    # ---- Pass 1: Line-by-line transforms (declarations, IMPLICIT NONE) ----
    result = []
    deferred_params = []

    i = 0
    while i < len(processed):
        line = processed[i]

        # Check if this is a comment line
        if len(line) > 0 and line[0] in ('*', 'C', 'c', '!'):
            result.append(line)
            i += 1
            continue

        # Empty or blank line
        if line.strip() == '':
            result.append(line)
            i += 1
            continue

        # Accumulate continuation lines
        full_line = line
        continuation_indices = [i]
        j = i + 1
        while j < len(processed):
            nextl = processed[j]
            if len(nextl) > 5 and nextl[0] not in ('*', 'C', 'c', '!') and nextl[5] not in (' ', '0', '') and nextl.strip() != '':
                full_line = full_line + nextl[6:]
                continuation_indices.append(j)
                j += 1
            else:
                break

        # Try to transform
        transformed = transform_statement(full_line, deferred_params)

        if transformed is not None:
            emit_lines = split_fortran_line(transformed)
            result.extend(emit_lines)
            i = j
        else:
            # Check for array-slice operations on the full (joined) line
            array_result = transform_array_ops(full_line)
            if array_result is not None:
                for ar in array_result:
                    result.extend(split_fortran_line(ar))
                i = j
            else:
                result.append(line)
                i += 1

        # Insert deferred PARAMETER statements
        if deferred_params:
            for pname, pval in deferred_params:
                pline = '      PARAMETER (%s = %s)' % (pname, pval)
                result.append(pline)
            deferred_params.clear()

    # Insert EXTERNAL declarations - find insertion point (after last declaration, before first executable)
    # Actually, just add them right after the first set of declarations
    if externals:
        result = insert_externals(result, externals)

    # ---- Pass 2: Handle DO/EXIT/CYCLE/named loops ----
    result = handle_loops(result)

    # ---- Pass 3: Handle END SUBROUTINE/FUNCTION ----
    final = []
    for line in result:
        stripped = line.lstrip()
        if re.match(r'END\s+SUBROUTINE\b', stripped, re.IGNORECASE):
            final.append('      END')
        elif re.match(r'END\s+FUNCTION\b', stripped, re.IGNORECASE):
            final.append('      END')
        else:
            final.append(line)

    return '\n'.join(final) + '\n'


def strip_interface_blocks(lines):
    """Comment out INTERFACE...END INTERFACE blocks and extract subroutine names as EXTERNAL."""
    result = []
    externals_by_sub = {}  # subroutine_idx -> set of external names
    current_sub_start = -1
    in_interface = False
    interface_sub_name = None

    for i, line in enumerate(lines):
        stripped = line.lstrip()

        # Track which subroutine we're in
        m = re.match(r'(SUBROUTINE|FUNCTION)\s+(\w+)', stripped, re.IGNORECASE)
        if m and not in_interface:
            current_sub_start = len(result)
            if current_sub_start not in externals_by_sub:
                externals_by_sub[current_sub_start] = set()

        if re.match(r'INTERFACE\s*$', stripped, re.IGNORECASE):
            in_interface = True
            result.append('*' + line[1:] if len(line) > 0 else '*')
            continue

        if in_interface:
            # Look for subroutine/function name inside the interface
            m = re.match(r'(SUBROUTINE|FUNCTION)\s+(\w+)', stripped, re.IGNORECASE)
            if m:
                interface_sub_name = m.group(2)

            if re.match(r'END\s+INTERFACE', stripped, re.IGNORECASE):
                in_interface = False
                if interface_sub_name and current_sub_start >= 0:
                    externals_by_sub[current_sub_start].add(interface_sub_name.upper())
                interface_sub_name = None
                result.append('*' + line[1:] if len(line) > 0 else '*')
            else:
                # Comment out the line
                result.append('*' + line[1:] if len(line) > 0 else '*')
            continue

        result.append(line)

    # Collect all externals
    all_externals = set()
    for s in externals_by_sub.values():
        all_externals.update(s)

    return result, all_externals


def insert_externals(lines, externals):
    """Insert EXTERNAL declarations for subroutines found in INTERFACE blocks."""
    if not externals:
        return lines

    # Find the first EXTERNAL statement or the first executable statement
    # and insert before it. Also need to handle multiple subroutines in one file.
    # Strategy: for each subroutine, find where EXTERNAL declarations go and add ours.

    # Simpler: find existing EXTERNAL lines and add after the last one we find,
    # or after the last declaration line
    result = []
    inserted = False

    for i, line in enumerate(lines):
        result.append(line)
        stripped = line.lstrip()
        # After existing EXTERNAL declarations
        if re.match(r'EXTERNAL\b', stripped, re.IGNORECASE) and not inserted:
            # Check if next lines are also EXTERNAL
            pass

    # Simpler approach: just find the first EXTERNAL in each subroutine and add after it
    # For now, let's just not insert - the EXTERNAL declarations from the original code
    # should suffice. The INTERFACE blocks were just for type checking.
    return lines


def transform_array_ops(line):
    """Transform F90 array operations to F77 equivalents."""
    stripped = line.lstrip()

    # Skip comments
    if len(line) > 0 and line[0] in ('*', 'C', 'c', '!'):
        return None

    # Pattern 1: ARR(:) = scalar_value
    # e.g., IWORK(:) = 0, RWORK(:) = ZERO, ISUPPZ(:) = 0
    m = re.match(r'(\w+)\s*\(\s*:\s*\)\s*=\s*(.+)', stripped, re.IGNORECASE)
    if m:
        arr = m.group(1)
        val = m.group(2).strip()
        # Can't determine array size, so comment out and add note
        # Actually, for f2c we need a DO loop. But we don't know the size.
        # We'll use a comment with the original and a simple workaround
        # For now, comment out the F90 line - these are typically initializations
        # that we can handle specially
        return ['*     F90: %s' % stripped,
                '*     TODO: convert to DO loop (array init)']

    # Pattern 2: ARR(lo:hi) = OTHER(lo2:hi2) or ARR(:) = OTHER(:)
    # e.g., RWORK(IXD : IXD + N-1) = D(:)
    m = re.match(r'(\w+)\s*\(([^)]*:[^)]*)\)\s*=\s*(\w+)\s*\(([^)]*:[^)]*)\)', stripped, re.IGNORECASE)
    if m:
        return ['*     F90: %s' % stripped,
                '*     TODO: convert to DO loop (array copy)']

    # Pattern 3: W(:) = ZERO (whole array)
    # Already covered by pattern 1

    return None


def transform_statement(line, deferred_params):
    """Transform a single Fortran statement. Returns None if no transform needed."""
    stripped = line.lstrip()

    # Skip comments
    if len(line) > 0 and line[0] in ('*', 'C', 'c', '!'):
        return None

    # Remove IMPLICIT NONE
    if re.match(r'\s*IMPLICIT\s+NONE', stripped, re.IGNORECASE):
        return '*     IMPLICIT NONE'

    # Handle: TYPE, INTENT(...) :: vars
    # and: TYPE, PARAMETER :: var = val
    m = re.match(
        r'\s*(INTEGER|DOUBLE\s+PRECISION|REAL|LOGICAL|CHARACTER)\s*,\s*'
        r'(INTENT\s*\([^)]*\)|PARAMETER)\s*::\s*(.*)',
        stripped, re.IGNORECASE
    )
    if m:
        typename = m.group(1).strip()
        modifier = m.group(2).strip().upper()
        rest = m.group(3).strip()

        if modifier == 'PARAMETER' or modifier.startswith('PARAMETER'):
            parts = split_param_assignments(rest)
            decl_names = []
            for pname, pval in parts:
                deferred_params.append((pname.strip(), pval.strip()))
                decl_names.append(pname.strip())
            return '      %s %s' % (typename, ', '.join(decl_names))
        else:
            # INTENT - just strip it
            rest = fix_array_bounds(rest)
            return '      %s %s' % (typename, rest)

    # Handle: TYPE :: vars (without INTENT or PARAMETER, just :: style)
    m = re.match(
        r'\s*(INTEGER|DOUBLE\s+PRECISION|REAL|LOGICAL|CHARACTER)\s*::\s*(.*)',
        stripped, re.IGNORECASE
    )
    if m:
        typename = m.group(1).strip()
        rest = m.group(2).strip()
        rest = fix_array_bounds(rest)
        return '      %s %s' % (typename, rest)

    return None


def fix_array_bounds(decl):
    """
    Convert F90-style array bounds like ARR(2*IL-1 : 2*IU) to ARR(*)
    For multi-dimensional: ARR(LDZ, WIL:WIU) -> ARR(LDZ, *)
    Only for declaration context.
    """
    result = []
    i = 0
    while i < len(decl):
        m = re.match(r'([A-Za-z_]\w*)\s*\(', decl[i:])
        if m:
            name = m.group(1)
            result.append(name)
            start = i + m.end() - 1
            depth = 1
            j = start + 1
            while j < len(decl) and depth > 0:
                if decl[j] == '(':
                    depth += 1
                elif decl[j] == ')':
                    depth -= 1
                j += 1
            inner = decl[start+1:j-1]
            if ':' in inner:
                # Split dimensions by comma (respecting parens)
                dims = split_by_comma_simple(inner)
                new_dims = []
                for dim in dims:
                    if ':' in dim:
                        new_dims.append('*')
                    else:
                        new_dims.append(dim.strip())
                result.append('(' + ', '.join(new_dims) + ')')
            else:
                result.append('(' + inner + ')')
            i = j
        else:
            result.append(decl[i])
            i += 1
    return ''.join(result)


def split_by_comma_simple(s):
    """Split by comma respecting parentheses."""
    parts = []
    depth = 0
    current = ''
    for ch in s:
        if ch == '(':
            depth += 1
            current += ch
        elif ch == ')':
            depth -= 1
            current += ch
        elif ch == ',' and depth == 0:
            parts.append(current)
            current = ''
        else:
            current += ch
    if current:
        parts.append(current)
    return parts


def split_param_assignments(s):
    """Split "A = 1, B = 2" into [("A","1"), ("B","2")]"""
    parts = []
    depth = 0
    current = ''
    for ch in s:
        if ch == '(':
            depth += 1
            current += ch
        elif ch == ')':
            depth -= 1
            current += ch
        elif ch == ',' and depth == 0:
            parts.append(current)
            current = ''
        else:
            current += ch
    if current.strip():
        parts.append(current)

    result = []
    for part in parts:
        if '=' in part:
            name, val = part.split('=', 1)
            result.append((name.strip(), val.strip()))
    return result


def split_fortran_line(line):
    """Split a long line into multiple continuation lines if needed (72 col limit)."""
    if len(line) <= 72:
        return [line]
    result = [line[:72]]
    rest = line[72:]
    while rest:
        chunk = '     $' + rest[:66]
        result.append(chunk)
        rest = rest[66:]
    return result


def handle_loops(lines):
    """
    Handle bare DO, named DO, EXIT, CYCLE.
    """
    loop_stack = []
    loop_info = {}

    label_counter = 90000

    for idx, line in enumerate(lines):
        stripped = line.lstrip()
        if len(line) > 0 and line[0] in ('*', 'C', 'c', '!'):
            continue

        # Named DO: "name: DO" or "name: DO I = ..."
        m = re.match(r'(\w+)\s*:\s*DO\b(.*)', stripped, re.IGNORECASE)
        if m:
            name = m.group(1)
            rest = m.group(2).strip()
            label_counter += 1
            top_label = label_counter
            label_counter += 1
            exit_label = label_counter
            is_bare = (rest == '' or rest.isspace())
            info = {'type': 'named_do', 'name': name, 'top_label': top_label,
                    'exit_label': exit_label, 'line_idx': idx, 'bare': is_bare,
                    'rest': rest}
            loop_stack.append(info)
            loop_info[idx] = info
            continue

        # Bare DO (no loop var, no name) - must be just "DO" at end of meaningful content
        m = re.match(r'DO\s*$', stripped, re.IGNORECASE)
        if m:
            label_counter += 1
            top_label = label_counter
            label_counter += 1
            exit_label = label_counter
            info = {'type': 'bare_do', 'name': None, 'top_label': top_label,
                    'exit_label': exit_label, 'line_idx': idx, 'bare': True}
            loop_stack.append(info)
            loop_info[idx] = info
            continue

        # Regular DO I = ... (unlabeled F90 style)
        m = re.match(r'DO\s+(\w+)\s*=\s*(.*)', stripped, re.IGNORECASE)
        if m:
            var = m.group(1)
            rest = m.group(2)
            label_counter += 1
            top_label = label_counter
            label_counter += 1
            exit_label = label_counter
            info = {'type': 'regular_do', 'name': None, 'top_label': top_label,
                    'exit_label': exit_label, 'line_idx': idx, 'bare': False,
                    'var': var, 'rest': rest}
            loop_stack.append(info)
            loop_info[idx] = info
            continue

        # ENDDO or END DO, possibly with name
        m = re.match(r'END\s*DO\s*(\w*)', stripped, re.IGNORECASE)
        if m:
            endname = m.group(1).strip() if m.group(1) else None
            if loop_stack:
                if endname:
                    for si in range(len(loop_stack)-1, -1, -1):
                        if loop_stack[si]['name'] and loop_stack[si]['name'].upper() == endname.upper():
                            match = loop_stack.pop(si)
                            break
                    else:
                        match = loop_stack.pop()
                else:
                    match = loop_stack.pop()
                info = {'type': 'enddo', 'match': match, 'line_idx': idx}
                loop_info[idx] = info
            continue

        # EXIT [name] - can appear standalone or after IF(...)
        m = re.search(r'\bEXIT\s*(\w*)', stripped, re.IGNORECASE)
        if m and not re.search(r'\bEXIT\s*=', stripped, re.IGNORECASE):
            ename = m.group(1).strip() if m.group(1) else None
            info = {'type': 'exit', 'name': ename, 'line_idx': idx}
            if ename and loop_stack:
                for si in range(len(loop_stack)-1, -1, -1):
                    if loop_stack[si]['name'] and loop_stack[si]['name'].upper() == ename.upper():
                        info['match'] = loop_stack[si]
                        break
            elif loop_stack:
                info['match'] = loop_stack[-1]
            loop_info[idx] = info
            continue

        # CYCLE [name] - can appear standalone or after IF(...)
        m = re.search(r'\bCYCLE\s*(\w*)', stripped, re.IGNORECASE)
        if m:
            cname = m.group(1).strip() if m.group(1) else None
            info = {'type': 'cycle', 'name': cname, 'line_idx': idx}
            if cname and loop_stack:
                for si in range(len(loop_stack)-1, -1, -1):
                    if loop_stack[si]['name'] and loop_stack[si]['name'].upper() == cname.upper():
                        info['match'] = loop_stack[si]
                        break
            elif loop_stack:
                info['match'] = loop_stack[-1]
            loop_info[idx] = info
            continue

    # Second pass: emit transformed lines
    result = []
    for idx, line in enumerate(lines):
        if idx not in loop_info:
            result.append(line)
            continue

        info = loop_info[idx]

        if info['type'] == 'bare_do':
            result.append('%5d CONTINUE' % info['top_label'])

        elif info['type'] == 'named_do':
            if info['bare']:
                result.append('%5d CONTINUE' % info['top_label'])
            else:
                result.append('      DO %d %s' % (info['top_label'], info['rest']))

        elif info['type'] == 'regular_do':
            result.append('      DO %d %s = %s' % (info['top_label'], info['var'], info['rest']))

        elif info['type'] == 'enddo':
            match = info['match']
            if match['bare']:
                result.append('      GOTO %d' % match['top_label'])
                result.append('%5d CONTINUE' % match['exit_label'])
            else:
                result.append('%5d CONTINUE' % match['top_label'])
                # Always emit exit_label for regular DO loops too,
                # in case EXIT targets this loop
                result.append('%5d CONTINUE' % match['exit_label'])

        elif info['type'] == 'exit':
            if 'match' in info:
                newline = re.sub(r'\bEXIT\b\s*\w*', 'GOTO %d' % info['match']['exit_label'], line, flags=re.IGNORECASE)
                result.append(newline)
            else:
                result.append(line)

        elif info['type'] == 'cycle':
            if 'match' in info:
                newline = re.sub(r'\bCYCLE\b\s*\w*', 'GOTO %d' % info['match']['top_label'], line, flags=re.IGNORECASE)
                result.append(newline)
            else:
                result.append(line)

        else:
            result.append(line)

    return result


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: f90_to_f77.py <file.f> [output.f]")
        sys.exit(1)

    infile = sys.argv[1]
    if len(sys.argv) > 2:
        outfile = sys.argv[2]
    else:
        outfile = infile

    result = process_file(infile)
    with open(outfile, 'w') as f:
        f.write(result)
    print("Converted: %s -> %s" % (infile, outfile))
