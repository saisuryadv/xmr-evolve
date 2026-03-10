#!/usr/bin/env python3
"""
Convert F90 array-slice operations to F77 DO loops.
Handles patterns like:
  ARR(lo:hi) = value
  ARR(lo:hi) = OTHER(lo2:hi2)
  ARR(lo:hi:stride) = value
  ARR(:) = value
  Z(1:N, WI:WJ) = ZERO  (2D)
"""

import re
import sys


def parse_slice(s):
    """Parse '1:N' or 'lo:hi:stride' or ':' into (lo, hi, stride)"""
    s = s.strip()
    parts = s.split(':')
    if len(parts) == 1:
        return None  # Not a slice, just an index
    lo = parts[0].strip() if parts[0].strip() else None
    hi = parts[1].strip() if len(parts) > 1 and parts[1].strip() else None
    stride = parts[2].strip() if len(parts) > 2 and parts[2].strip() else None
    return (lo, hi, stride)


def is_slice_expr(s):
    """Check if expression contains : indicating a slice"""
    # Careful: could be inside function calls
    depth = 0
    for ch in s:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        elif ch == ':' and depth == 0:
            return True
    return False


def find_matching_paren(s, start):
    """Find matching closing paren starting from opening paren at start"""
    depth = 1
    i = start + 1
    while i < len(s) and depth > 0:
        if s[i] == '(':
            depth += 1
        elif s[i] == ')':
            depth -= 1
        i += 1
    return i - 1  # position of closing paren


def process_line(line, linenum, filename):
    """Process a single line. Returns list of replacement lines, or None if no change."""
    stripped = line.strip()

    # Skip comments
    if not stripped or stripped[0] in ('*', 'C', 'c', '!'):
        return None

    # Skip lines with :: (declarations), INTENT, INTERFACE
    upper = stripped.upper()
    if '::' in stripped or 'INTENT' in upper or 'INTERFACE' in upper:
        return None

    # Check for array slice assignment: LHS(...:...) = RHS
    # Need to find = sign that's the assignment (not inside parens)
    eq_pos = find_assignment_eq(stripped)
    if eq_pos is None:
        return None

    lhs = stripped[:eq_pos].strip()
    rhs = stripped[eq_pos+1:].strip()

    # Check if LHS has array slice
    if not has_slice(lhs):
        # Check if RHS has array slice (but LHS doesn't)
        # e.g., RWORK(IXROOT : IXROOT+RSTRIDER-1) = ROOTR
        # Wait, that's LHS with slice. Let me re-examine.
        return None

    # Parse LHS
    # Could be: ARR(slice) or ARR(slice, slice) or ARR(idx, slice)
    m = re.match(r'(\w+)\s*\((.+)\)', lhs)
    if not m:
        return None

    arrname = m.group(1)
    args_str = m.group(2)

    # Split args by comma (respecting parens)
    args = split_by_comma(args_str)

    # Determine dimensionality and which dims are slices
    dim_info = []
    for arg in args:
        sl = parse_slice(arg)
        if sl is not None:
            dim_info.append(('slice', sl))
        else:
            dim_info.append(('index', arg.strip()))

    # Parse RHS similarly
    rhs_m = re.match(r'(\w+)\s*\((.+)\)', rhs)
    rhs_is_array = False
    rhs_dims = []
    if rhs_m and has_slice(rhs):
        rhs_arrname = rhs_m.group(1)
        rhs_args = split_by_comma(rhs_m.group(2))
        for arg in rhs_args:
            sl = parse_slice(arg)
            if sl is not None:
                rhs_dims.append(('slice', sl))
            else:
                rhs_dims.append(('index', arg.strip()))
        rhs_is_array = True

    # Get indentation
    indent = len(line) - len(line.lstrip())
    pad = ' ' * indent

    # Generate DO loop(s)
    # Use loop variable names ITMP1, ITMP2, etc.
    loop_vars = []
    slice_dims = [(i, info) for i, info in enumerate(dim_info) if info[0] == 'slice']

    if len(slice_dims) == 0:
        return None

    lines_out = []
    loop_var_names = ['IXF77A', 'IXF77B', 'IXF77C']

    # Build nested DO loops (outermost first)
    for li, (di, (_, (lo, hi, stride))) in enumerate(slice_dims):
        lv = loop_var_names[li]
        loop_vars.append((di, lv, lo, hi, stride))

        if lo is None:
            lo = '1'
        if hi is None:
            # We don't know the size - use a placeholder
            hi = '???'

        if stride:
            lines_out.append('%sDO %s = %s, %s, %s' % (pad, lv, lo, hi, stride))
        else:
            lines_out.append('%sDO %s = %s, %s' % (pad, lv, lo, hi))

    # Build the assignment
    # LHS index expression
    lhs_indices = []
    lv_idx = 0
    for i, info in enumerate(dim_info):
        if info[0] == 'slice':
            lhs_indices.append(loop_vars[lv_idx][1])  # loop var name
            lv_idx += 1
        else:
            lhs_indices.append(info[1])

    lhs_expr = '%s(%s)' % (arrname, ', '.join(lhs_indices))

    # RHS expression
    if rhs_is_array:
        rhs_indices = []
        lv_idx = 0
        for i, info in enumerate(rhs_dims):
            if info[0] == 'slice':
                _, (rlo, rhi, rstride) = info
                # Map LHS loop var to RHS range
                # If LHS goes lo..hi and RHS goes rlo..rhi, then
                # RHS index = rlo + (lv - lo) if lo is not None
                lhs_di, lhs_lv, lhs_lo, lhs_hi, lhs_stride = loop_vars[lv_idx]
                if rlo and lhs_lo:
                    if rlo == lhs_lo:
                        rhs_indices.append(lhs_lv)
                    else:
                        rhs_indices.append('%s + (%s) - (%s)' % (rlo, lhs_lv, lhs_lo if lhs_lo else '1'))
                elif rlo:
                    rhs_indices.append('%s + %s - 1' % (rlo, lhs_lv))
                else:
                    rhs_indices.append(lhs_lv)
                lv_idx += 1
            else:
                rhs_indices.append(info[1])

        rhs_expr = '%s(%s)' % (rhs_arrname, ', '.join(rhs_indices))
    else:
        rhs_expr = rhs

    lines_out.append('%s   %s = %s' % (pad, lhs_expr, rhs_expr))

    # Close loops (innermost first)
    for li in range(len(slice_dims)-1, -1, -1):
        lines_out.append('%sENDDO' % pad)

    return lines_out


def find_assignment_eq(s):
    """Find the position of = sign that is the assignment operator (not inside parens, not ==, not <=, etc.)"""
    depth = 0
    i = 0
    while i < len(s):
        if s[i] == '(':
            depth += 1
        elif s[i] == ')':
            depth -= 1
        elif s[i] == '=' and depth == 0:
            # Check it's not == or <= or >= or /=
            if i > 0 and s[i-1] in ('<', '>', '/', '='):
                i += 1
                continue
            if i + 1 < len(s) and s[i+1] == '=':
                i += 2
                continue
            return i
        i += 1
    return None


def has_slice(s):
    """Check if string contains array slice notation (: inside parens)"""
    depth = 0
    in_parens = False
    for ch in s:
        if ch == '(':
            depth += 1
            in_parens = True
        elif ch == ')':
            depth -= 1
            if depth == 0:
                in_parens = False
        elif ch == ':' and in_parens and depth > 0:
            return True
    return False


def split_by_comma(s):
    """Split string by commas, respecting parentheses"""
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
    return parts


def process_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    output = []
    changes = 0
    needs_itmp = False

    for i, line in enumerate(lines):
        result = process_line(line.rstrip('\n'), i+1, filepath)
        if result is not None:
            # Check for ??? placeholders
            has_unknown = any('???' in r for r in result)
            if has_unknown:
                output.append(line.rstrip('\n'))
                print("  WARNING line %d: unknown array bounds: %s" % (i+1, line.strip()))
            else:
                output.extend(result)
                changes += 1
                needs_itmp = True
        else:
            output.append(line.rstrip('\n'))

    if changes > 0:
        # Add IXF77A etc. to local declarations
        # Find the last INTEGER declaration line and add after it
        insert_idx = None
        for j in range(len(output)):
            stripped = output[j].lstrip()
            if stripped.upper().startswith('INTEGER') and '::' not in stripped and 'INTENT' not in stripped.upper():
                insert_idx = j

        if insert_idx is not None and needs_itmp:
            output.insert(insert_idx + 1, '      INTEGER IXF77A, IXF77B, IXF77C')

        with open(filepath, 'w') as f:
            f.write('\n'.join(output) + '\n')
        print("  Fixed %d array-slice operations" % changes)
    else:
        print("  No array-slice operations found")

    return changes


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: fix_array_slices.py <file.f> [file2.f ...]")
        sys.exit(1)

    for filepath in sys.argv[1:]:
        print("Processing: %s" % filepath)
        process_file(filepath)
