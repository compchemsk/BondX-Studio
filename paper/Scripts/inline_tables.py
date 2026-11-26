#!/usr/bin/env python3
"""inline_tables.py
Writes ESI_V3.tex from ESI_V3_src.tex. Every line \\input{tables/NAME.tex} is replaced by the text of tables/NAME.tex.
The tables come from scripts 11 and 26. The result holds all tables inside one .tex file.
Usage: python3 inline_tables.py /path/to/paper_folder
"""
import os, re, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
SRC, OUT = os.path.join(ROOT, 'ESI_V3_src.tex'), os.path.join(ROOT, 'ESI_V3.tex')
pat = re.compile(r'^\\input\{(tables/[A-Za-z0-9_]+)\.tex\}\s*$')
out, n = [], 0
for line in open(SRC, encoding='utf-8'):
    m = pat.match(line.strip())
    if not m:
        out.append(line)
        continue
    path = os.path.join(ROOT, m.group(1) + '.tex')
    body = open(path, encoding='utf-8').read().rstrip('\n')
    out.append('%% ---- table from %s.tex (generated) ----\n%s\n%% ---- end of %s.tex ----\n' % (m.group(1), body, m.group(1)))
    n += 1
text = ''.join(out)
text = text.replace('%% ESI_V3_src.tex -- source of ESI_V3.tex', '%% ESI_V3.tex -- written by Scripts/inline_tables.py from ESI_V3_src.tex', 1)
assert not any(l.lstrip().startswith('\\input{tables/') for l in text.splitlines()), 'a table input was not replaced'
open(OUT, 'w', encoding='utf-8').write(text)
print('inlined %d tables into %s' % (n, OUT))
