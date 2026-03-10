#!/usr/bin/env python3
"""
Preprocess Fortran 90 files to Fortran 77 compatible format for f2c.
Handles the specific F90 features used in XMR and modern LAPACK files:
- INTENT(IN/OUT/INOUT) declarations
- IMPLICIT NONE
- TYPE, PARAMETER :: NAME = VALUE
- END DO / END IF / END SUBROUTINE / END FUNCTION
- DO ... END DO (without labels)
- Declarations mixed with PARAMETER (reorder to F77 order)
"""

import re
import sys
import os


def preprocess_f90_to_f77(lines):
    """Convert F90 lines to F77 compatible lines.

    Strategy: Two-pass approach.
    Pass 1: Collect all lines, converting F90 syntax to F77 equivalents.
    Pass 2: Reorder declarations to ensure all come before executables.
    """
    # Pass 1: syntax conversion
    converted = []
    do_stack = []
    next_label = 9000

    i = 0
    while i < len(lines):
        line = lines[i]

        # Preserve blank lines
        if not line.strip():
            converted.append(line)
            i += 1
            continue

        # Preserve comment lines (column 1 is C, c, *, or !)
        if len(line) > 0 and line[0] in ('C', 'c', '*', '!'):
            converted.append(line)
            i += 1
            continue

        # Get content after column 6
        if len(line) >= 6:
            col1_6 = line[:6]
            rest = line[6:] if len(line) > 6 else ''
        else:
            col1_6 = line.ljust(6)
            rest = ''

        # Check if continuation line (non-blank, non-zero in column 6)
        is_continuation = len(line) > 5 and line[5] not in (' ', '0', '\t') and col1_6.strip() == ''

        rest_stripped = rest.strip()
        rest_upper = rest_stripped.upper()

        # Skip IMPLICIT NONE
        if rest_upper == 'IMPLICIT NONE':
            converted.append('C     IMPLICIT NONE')
            i += 1
            continue

        # Handle TYPE, PARAMETER :: NAME = VALUE
        # e.g., INTEGER, PARAMETER :: TOLFAC = 3
        # → INTEGER TOLFAC
        #   PARAMETER (TOLFAC = 3)
        param_match = re.match(
            r'\s*(INTEGER|DOUBLE\s+PRECISION|REAL|LOGICAL|CHARACTER)(\*\d+)?\s*,\s*PARAMETER\s*::\s*(\w+)\s*=\s*(.*)',
            rest, re.IGNORECASE
        )
        if param_match:
            dtype = param_match.group(1).upper()
            size = param_match.group(2) or ''
            name = param_match.group(3)
            value = param_match.group(4).strip()
            converted.append(f'      {dtype}{size} {name}')
            converted.append(f'      PARAMETER ({name} = {value})')
            i += 1
            continue

        # Handle TYPE, INTENT(...) :: varlist
        intent_match = re.match(
            r'\s*(INTEGER|DOUBLE\s+PRECISION|REAL|LOGICAL|CHARACTER)(\*\d+)?\s*,\s*INTENT\s*\(\s*\w+\s*\)\s*::\s*(.*)',
            rest, re.IGNORECASE
        )
        if intent_match:
            dtype = intent_match.group(1).upper()
            size = intent_match.group(2) or ''
            varlist = intent_match.group(3).strip()
            new_line = f'      {dtype}{size} {varlist}'
            if len(new_line) > 72:
                parts = varlist.split(',')
                converted.append(f'      {dtype}{size} {parts[0].strip()}')
                for p in parts[1:]:
                    converted.append(f'     $   , {p.strip()}')
            else:
                converted.append(new_line)
            i += 1
            continue

        # Handle standalone INTENT declarations
        if re.match(r'\s*INTENT\s*\(\s*\w+\s*\)\s*::', rest, re.IGNORECASE):
            converted.append('C' + line[1:])
            i += 1
            continue

        # Handle END DO → CONTINUE with label
        if rest_upper in ('END DO', 'ENDDO'):
            if do_stack:
                label = do_stack.pop()
                converted.append(f' {label} CONTINUE')
            else:
                converted.append('C     END DO (unmatched)')
            i += 1
            continue

        # Handle END IF
        if rest_upper in ('END IF', 'ENDIF'):
            converted.append('      ENDIF')
            i += 1
            continue

        # Handle END SUBROUTINE / END FUNCTION
        if re.match(r'\s*END\s+(SUBROUTINE|FUNCTION)', rest, re.IGNORECASE):
            converted.append('      END')
            i += 1
            continue

        # Handle unlabeled DO var = ...
        do_match = re.match(r'(\s*)DO\s+([A-Za-z]\w*)\s*=\s*(.*)', rest, re.IGNORECASE)
        if do_match and not is_continuation:
            label_in_col = col1_6.strip()
            if not label_in_col:
                # Check if DO <number> form (already labeled)
                do_num = re.match(r'(\s*)DO\s+(\d+)\s+', rest, re.IGNORECASE)
                if do_num:
                    converted.append(line)
                    i += 1
                    continue
                # Unlabeled DO
                label = next_label
                next_label += 1
                do_stack.append(label)
                var = do_match.group(2)
                bounds = do_match.group(3)
                converted.append(f'      DO {label} {var} = {bounds}')
                i += 1
                continue

        # Handle DO WHILE
        if re.match(r'\s*DO\s+WHILE', rest, re.IGNORECASE) and not is_continuation:
            label = next_label
            next_label += 1
            do_stack.append(label)
            # Rewrite with label
            while_part = re.sub(r'^DO\s+WHILE', f'DO {label} WHILE', rest.strip(), flags=re.IGNORECASE)
            converted.append(f'      {while_part}')
            i += 1
            continue

        # Default: pass through
        converted.append(line)
        i += 1

    # Pass 2: Reorder declarations before executables within each program unit
    result = reorder_declarations(converted)
    return result


def reorder_declarations(lines):
    """Ensure all declarations come before executable statements.

    F90 allows mixing, F77 does not. We identify program unit boundaries
    and within each, move declarations before executables.
    """
    # For simplicity and safety, we'll just ensure PARAMETER statements
    # come right after their type declarations, before any executable code.
    # The main issue is that our converted PARAMETER lines might be
    # interleaved with other declarations.

    result = []
    in_decl_section = True
    decl_lines = []
    exec_lines = []
    unit_started = False

    for line in lines:
        stripped = line.strip().upper()

        # Detect program unit start
        if re.match(r'\s*SUBROUTINE\s+', stripped) or re.match(r'\s*FUNCTION\s+', stripped):
            # Flush previous unit
            if unit_started:
                result.extend(decl_lines)
                result.extend(exec_lines)
            decl_lines = [line]
            exec_lines = []
            in_decl_section = True
            unit_started = True
            continue

        if not unit_started:
            result.append(line)
            continue

        # Comments and blank lines go to current section
        if not stripped or (len(line) > 0 and line[0] in ('C', 'c', '*', '!')):
            if in_decl_section:
                decl_lines.append(line)
            else:
                exec_lines.append(line)
            continue

        # Check if this is a declaration line
        is_decl = False
        for pat in [
            r'INTEGER', r'DOUBLE\s+PRECISION', r'REAL', r'LOGICAL', r'CHARACTER',
            r'PARAMETER\s*\(', r'EXTERNAL', r'INTRINSIC', r'DIMENSION',
            r'COMMON', r'EQUIVALENCE', r'DATA', r'SAVE',
            r'C\s+IMPLICIT',  # Our commented-out IMPLICIT NONE
        ]:
            if re.match(r'\s*' + pat, stripped):
                is_decl = True
                break

        # Continuation lines belong to whatever section we're in
        if len(line) > 5 and line[5] not in (' ', '0', '\t') and line[:5].strip() == '':
            if in_decl_section:
                decl_lines.append(line)
            else:
                exec_lines.append(line)
            continue

        if is_decl:
            if in_decl_section:
                decl_lines.append(line)
            else:
                # Declaration after executable — move it to decl section
                decl_lines.append(line)
            continue

        # END statement
        if stripped == 'END':
            if in_decl_section:
                result.extend(decl_lines)
                decl_lines = []
            else:
                result.extend(decl_lines)
                result.extend(exec_lines)
                decl_lines = []
                exec_lines = []
            result.append(line)
            in_decl_section = True
            unit_started = False
            continue

        # Executable statement
        if in_decl_section:
            in_decl_section = False
        exec_lines.append(line)

    # Flush last unit
    if unit_started:
        result.extend(decl_lines)
        result.extend(exec_lines)

    return result


def convert_file(input_path, output_path):
    """Convert a single F90 file to F77."""
    with open(input_path, 'r') as f:
        lines = f.readlines()
    lines = [l.rstrip('\n').rstrip('\r') for l in lines]
    result = preprocess_f90_to_f77(lines)
    with open(output_path, 'w') as f:
        for line in result:
            f.write(line + '\n')
    return len(result)


def main():
    if len(sys.argv) < 2:
        print("Usage: f90_to_f77.py <input.f> [output.f]")
        print("       f90_to_f77.py <directory> [output_dir]")
        sys.exit(1)

    path = sys.argv[1]

    if os.path.isdir(path):
        out_dir = sys.argv[2] if len(sys.argv) > 2 else path + '_f77'
        os.makedirs(out_dir, exist_ok=True)
        for fname in sorted(os.listdir(path)):
            if fname.endswith('.f'):
                inp = os.path.join(path, fname)
                outp = os.path.join(out_dir, fname)
                n = convert_file(inp, outp)
                print(f"  {fname}: {n} lines")
        print(f"Output in: {out_dir}")
    else:
        out_path = sys.argv[2] if len(sys.argv) > 2 else path.replace('.f', '_f77.f')
        n = convert_file(path, out_path)
        print(f"Converted {path} → {out_path} ({n} lines)")


if __name__ == '__main__':
    main()
