#!/usr/bin/env python3
"""
Extract all distinct versions of bidiag_svd.h from the JSONL transcript.
Version 3: Properly handles failed edits and file state tracking.

Strategy:
- Parse all events (Read, Write, Edit, Eval) in order
- Track file state from Read results (which show actual disk state)
- For bidiag_svd.h: start from Read result, apply edits, snapshot after evals
- For dlaxre.c/dlaxre.f: Read results after edits show the modified state
- Mark edit as failed only if tool_result says so (line 375 is the only failure)
"""
import json
import re
import os
import sys

JSONL_PATH = '/Users/saisurya/.claude/projects/-Users-saisurya-MRRR-xmr-evolve/b05f34bb-76b2-42f5-8184-65ab5d3937ea.jsonl'
OUTPUT_DIR = '/Users/saisurya/MRRR/xmr-evolve/recovered'
OTHER_DIR = os.path.join(OUTPUT_DIR, 'other_files')
BIDIAG_PATH = '/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_035055/agent_79337_0001/program/bidiag_svd.h'
DLAXRE_C_PATH = '/Users/saisurya/MRRR/xmr-evolve/src/xmr_c/dlaxre.c'
DLAXRE_F_PATH = '/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_035055/agent_79337_0001/xmr_src/dlaxre.f'

os.makedirs(OTHER_DIR, exist_ok=True)

# Load all lines
with open(JSONL_PATH) as f:
    all_lines = [json.loads(l) for l in f]

print(f"Total JSONL lines: {len(all_lines)}")

# Known failed edits (tool_result showed error)
FAILED_EDITS = {375}  # line 375: "File has not been read yet"

def strip_read_output(text):
    """Strip line number prefixes from Read tool output.
    Format: '     1→content' where → is U+2192."""
    raw_lines = text.split('\n')
    clean = []
    for l in raw_lines:
        arrow_idx = l.find('\u2192')
        if arrow_idx >= 0:
            prefix = l[:arrow_idx].strip()
            if prefix.isdigit():
                clean.append(l[arrow_idx+1:])
            else:
                clean.append(l)
        else:
            clean.append(l)
    return '\n'.join(clean)


def get_tool_result_for_id(tool_use_id, start_line):
    """Find the tool_result for a given tool_use_id."""
    for i in range(start_line, min(start_line + 15, len(all_lines))):
        obj = all_lines[i]
        if obj.get('type') == 'user':
            msg = obj.get('message', {})
            content = msg.get('content', '')
            if isinstance(content, list):
                for block in content:
                    if isinstance(block, dict) and block.get('type') == 'tool_result':
                        if block.get('tool_use_id') == tool_use_id:
                            return block.get('content', '')
    return None


def find_read_result(path, after_line):
    """Find a Read tool call for path after given line, return stripped content."""
    for i in range(after_line, len(all_lines)):
        obj = all_lines[i]
        if obj.get('type') == 'assistant':
            msg = obj.get('message', {})
            content = msg.get('content', [])
            if isinstance(content, list):
                for block in content:
                    if isinstance(block, dict) and block.get('type') == 'tool_use' and block.get('name') == 'Read':
                        fp = block.get('input', {}).get('file_path', '')
                        if fp == path:
                            tool_id = block.get('id', '')
                            result = get_tool_result_for_id(tool_id, i+1)
                            if isinstance(result, str) and len(result) > 100 and '\u2192' in result:
                                return strip_read_output(result), i
    return None, None


# Collect ALL events in order
events = []

for i, obj in enumerate(all_lines):
    tp = obj.get('type', '?')

    if tp == 'assistant':
        msg = obj.get('message', {})
        content = msg.get('content', [])
        if isinstance(content, list):
            for block in content:
                if isinstance(block, dict) and block.get('type') == 'tool_use':
                    tool = block.get('name', '?')
                    inp = block.get('input', {})
                    tool_id = block.get('id', '')

                    if tool == 'Write':
                        events.append((i, 'write', {
                            'file_path': inp.get('file_path', ''),
                            'content': inp.get('content', ''),
                            'tool_id': tool_id,
                        }))
                    elif tool == 'Edit':
                        events.append((i, 'edit', {
                            'file_path': inp.get('file_path', ''),
                            'old_string': inp.get('old_string', ''),
                            'new_string': inp.get('new_string', ''),
                            'replace_all': inp.get('replace_all', False),
                            'tool_id': tool_id,
                        }))
                    elif tool == 'Read':
                        fp = inp.get('file_path', '')
                        events.append((i, 'read', {
                            'file_path': fp,
                            'tool_id': tool_id,
                        }))

    elif tp == 'user':
        msg = obj.get('message', {})
        content = msg.get('content', '')
        if isinstance(content, list):
            for block in content:
                if isinstance(block, dict) and block.get('type') == 'tool_result':
                    result_content = block.get('content', '')
                    result_texts = []
                    if isinstance(result_content, str):
                        result_texts = [result_content]
                    elif isinstance(result_content, list):
                        for rb in result_content:
                            if isinstance(rb, dict):
                                result_texts.append(rb.get('text', ''))

                    for text in result_texts:
                        if '/289' in text:
                            m = re.search(r'Pass:\s*(\d+)/289', text)
                            if m:
                                pass_rate = int(m.group(1))
                                score_m = re.search(r'composite_score=([\d.]+)', text)
                                score = score_m.group(1) if score_m else '?'
                                scaling_m = re.search(r'pass_worst_scaling=([\d.]+)', text)
                                scaling = scaling_m.group(1) if scaling_m else '?'
                                worst_m = re.search(r'pass_worst_scaling=[\d.]+ \((\S+)', text)
                                worst_matrix = worst_m.group(1) if worst_m else '?'
                                pr_m = re.search(r'pass_rate=([\d.]+)', text)
                                pr = pr_m.group(1) if pr_m else '?'
                                res_m = re.search(r'avg:\s*res=([\d.]+)', text)
                                res_avg = res_m.group(1) if res_m else '?'
                                orthoU_m = re.search(r'orthoU=([\d.]+)', text)
                                orthoU = orthoU_m.group(1) if orthoU_m else '?'
                                orthoV_m = re.search(r'orthoV=([\d.]+)', text)
                                orthoV = orthoV_m.group(1) if orthoV_m else '?'

                                events.append((i, 'eval', {
                                    'pass_rate': pass_rate,
                                    'score': score,
                                    'scaling': scaling,
                                    'worst_matrix': worst_matrix,
                                    'pass_rate_frac': pr,
                                    'res_avg': res_avg,
                                    'orthoU': orthoU,
                                    'orthoV': orthoV,
                                    'full_text': text,
                                }))

# Get the initial bidiag_svd.h from the first successful Read
original_bidiag, orig_line = find_read_result(BIDIAG_PATH, 0)
print(f"Original bidiag_svd.h: {len(original_bidiag) if original_bidiag else 0} chars from line {orig_line}")

# Get initial dlaxre.f
original_dlaxre_f, _ = find_read_result(DLAXRE_F_PATH, 0)
print(f"Original dlaxre.f: {len(original_dlaxre_f) if original_dlaxre_f else 0} chars")

# Get initial dlaxre.c
original_dlaxre_c, _ = find_read_result(DLAXRE_C_PATH, 0)
print(f"Original dlaxre.c: {len(original_dlaxre_c) if original_dlaxre_c else 0} chars")


def apply_edit(content, old_string, new_string, replace_all, line_num):
    """Apply edit to content. Returns (new_content, success, error_msg)."""
    if content is None:
        return None, False, "content is None"
    if old_string not in content:
        return content, False, f"old_string not found ({len(old_string)} chars)"
    if replace_all:
        return content.replace(old_string, new_string), True, None
    else:
        return content.replace(old_string, new_string, 1), True, None


# Process events
tracked = {}
tracked[BIDIAG_PATH] = original_bidiag
if original_dlaxre_f:
    tracked[DLAXRE_F_PATH] = original_dlaxre_f
if original_dlaxre_c:
    tracked[DLAXRE_C_PATH] = original_dlaxre_c

versions = []
version_num = 0
pending_changes = []
edit_failures = []

for idx, (line_num, etype, data) in enumerate(events):
    if etype == 'read':
        fp = data['file_path']
        # When a Read happens, it reflects the actual state on disk.
        # For dlaxre.c, which is edited on the actual disk, the Read after edits
        # shows the updated state. We use this to sync our tracking.
        # However, for bidiag_svd.h in the agent workspace, our tracking should
        # already be correct. Let's just use Reads as a reference check.
        pass

    elif etype == 'write':
        fp = data['file_path']
        content = data['content']
        tracked[fp] = content
        pending_changes.append(f"Line {line_num}: WRITE {os.path.basename(fp)} ({len(content)} chars)")
        if fp == BIDIAG_PATH:
            print(f"\nWRITE bidiag_svd.h at line {line_num} ({len(content)} chars)")

    elif etype == 'edit':
        fp = data['file_path']

        # Skip known failed edits
        if line_num in FAILED_EDITS:
            msg = f"Line {line_num}: EDIT {os.path.basename(fp)} SKIPPED (known failure)"
            edit_failures.append(msg)
            print(f"  SKIP edit at line {line_num} (known failure)")
            continue

        if fp not in tracked:
            # Try to find content from a preceding Read
            content, _ = find_read_result(fp, 0)
            if content:
                tracked[fp] = content
                print(f"  Auto-loaded {os.path.basename(fp)} from Read ({len(content)} chars)")
            else:
                msg = f"Line {line_num}: EDIT {fp} - file not tracked"
                edit_failures.append(msg)
                print(f"  WARNING: {msg}")
                continue

        new_content, success, err = apply_edit(
            tracked[fp], data['old_string'], data['new_string'],
            data['replace_all'], line_num
        )

        if success:
            tracked[fp] = new_content
            basename = os.path.basename(fp)
            old_len = len(data['old_string'])
            new_len = len(data['new_string'])
            pending_changes.append(f"Line {line_num}: EDIT {basename} (old={old_len}, new={new_len})")
        else:
            msg = f"Line {line_num}: EDIT {os.path.basename(fp)} FAILED: {err}"
            edit_failures.append(msg)
            pending_changes.append(msg)
            print(f"  EDIT FAIL at line {line_num}: {err}")

            # Since we know all edits succeeded on the actual agent (except line 375),
            # if our tracking fails, resync from the next Read result
            read_content, read_line = find_read_result(fp, line_num)
            if read_content:
                print(f"  Resyncing {os.path.basename(fp)} from Read at line {read_line}")
                tracked[fp] = read_content

    elif etype == 'eval':
        if BIDIAG_PATH in tracked and tracked[BIDIAG_PATH]:
            pr = data['pass_rate']
            fname = f'version_{version_num:03d}_pass{pr:03d}.h'
            version_path = os.path.join(OUTPUT_DIR, fname)
            with open(version_path, 'w') as f:
                f.write(tracked[BIDIAG_PATH])

            versions.append({
                'num': version_num,
                'pass_rate': pr,
                'score': data['score'],
                'scaling': data['scaling'],
                'worst_matrix': data.get('worst_matrix', '?'),
                'pass_rate_frac': data.get('pass_rate_frac', '?'),
                'res_avg': data.get('res_avg', '?'),
                'orthoU': data.get('orthoU', '?'),
                'orthoV': data.get('orthoV', '?'),
                'line_num': line_num,
                'content_len': len(tracked[BIDIAG_PATH]),
                'changes': list(pending_changes),
                'fname': fname,
            })

            print(f"  VERSION {version_num:03d} Pass:{pr}/289 score={data['score']} scaling={data['scaling']} ({len(tracked[BIDIAG_PATH])} chars)")
            version_num += 1
            pending_changes = []

# Save other files
print(f"\n=== Saving other files ===")
saved_others = []
for fp, content in tracked.items():
    if fp == BIDIAG_PATH:
        continue
    basename = os.path.basename(fp)
    if 'xmr_src' in fp or 'xmr_c' in fp:
        parent = os.path.basename(os.path.dirname(fp))
        out_name = f"{parent}_{basename}"
    else:
        out_name = basename
    out_path = os.path.join(OTHER_DIR, out_name)
    with open(out_path, 'w') as f:
        f.write(content)
    saved_others.append((out_name, len(content), fp))
    print(f"  {out_name} ({len(content)} chars)")

# Write SUMMARY.md
print(f"\n=== Writing SUMMARY.md ===")
L = []
L.append("# Recovered Versions of bidiag_svd.h")
L.append("")
L.append("## Source")
L.append("")
L.append("- **JSONL**: `b05f34bb-76b2-42f5-8184-65ab5d3937ea.jsonl`")
L.append("- **Agent workspace**: `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_035055/agent_79337_0001/`")
L.append(f"- **Total versions extracted**: {len(versions)}")
L.append(f"- **Pass rates seen**: {', '.join(str(v['pass_rate']) for v in versions)}")
L.append("")
L.append("## Overview")
L.append("")
L.append("The agent explored a C++ bidiagonal SVD implementation using the Willems XMR")
L.append("(improved MR3) eigensolver. It started from a template `bidiag_svd.h` that calls")
L.append("DSTEXR (from the XMR library) on TGK tridiagonal matrices, then extracts U/V")
L.append("from eigenvectors. Over 32 edits and 19 evaluations, it improved from 114/289 to")
L.append("a peak of 207/289 passing tests.")
L.append("")

L.append("## Version Table")
L.append("")
L.append("| # | File | Pass | Score | Scaling | Worst Matrix | Size |")
L.append("|---|------|------|-------|---------|--------------|------|")

for v in versions:
    L.append(f"| {v['num']:03d} | `{v['fname']}` | {v['pass_rate']}/289 | {v['score']} | {v['scaling']} | {v['worst_matrix']} | {v['content_len']} |")

L.append("")
L.append("## Version Details")
L.append("")

for v in versions:
    L.append(f"### Version {v['num']:03d} -- Pass: {v['pass_rate']}/289")
    L.append("")
    L.append(f"- **File**: `{v['fname']}`")
    L.append(f"- **Composite Score**: {v['score']}")
    L.append(f"- **Worst Scaling**: {v['scaling']} ({v['worst_matrix']})")
    L.append(f"- **Avg Metrics**: res={v['res_avg']}, orthoU={v['orthoU']}, orthoV={v['orthoV']}")
    L.append(f"- **JSONL eval line**: {v['line_num']}")
    L.append(f"- **Content size**: {v['content_len']} chars")
    if v['changes']:
        L.append(f"- **Changes**:")
        for c in v['changes']:
            L.append(f"  - {c}")
    L.append("")

L.append("## Other Recovered Files")
L.append("")
L.append("Located in `other_files/`:")
L.append("")
for name, size, orig_path in saved_others:
    L.append(f"- `{name}` ({size} chars) -- from `{orig_path}`")

if edit_failures:
    L.append("")
    L.append("## Edit Failures / Notes")
    L.append("")
    for msg in edit_failures:
        L.append(f"- {msg}")

L.append("")
L.append("## Key Observations")
L.append("")
L.append("1. **Best accuracy**: Version 017 achieved 207/289 (71.6%) but with terrible scaling (19.50x)")
L.append("2. **Best balanced**: Version 010 and 013 achieved 138/289 with good scaling (4.74x and 4.29x)")
L.append("3. **Regression in v018**: Dropping from 207 to 171 after the last edit suggests the final")
L.append("   change hurt some test cases")
L.append("4. **Score of 5.0000**: Appears when scaling exceeds threshold (penalty)")
L.append("5. The agent modified `dlaxre.c` (the C wrapper for DSTEXR) to change eigenvalue")
L.append("   extraction behavior, and `dlaxre.f` (Fortran source) for the same purpose")

summary_path = os.path.join(OUTPUT_DIR, 'SUMMARY.md')
with open(summary_path, 'w') as f:
    f.write('\n'.join(L))

print(f"\nSUMMARY.md written to {summary_path}")
print(f"\nFinal summary:")
print(f"  Versions: {len(versions)}")
print(f"  Other files: {len(saved_others)}")
print(f"  Edit failures: {len(edit_failures)}")
for v in versions:
    print(f"  v{v['num']:03d}: Pass {v['pass_rate']}/289, score={v['score']}, scaling={v['scaling']}")
