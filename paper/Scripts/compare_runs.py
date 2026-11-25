#!/usr/bin/env python3
"""compare_runs.py  RUN_A RUN_B
Compares every CSV and JSON file of results/ between two pipeline runs.
Timing columns (ms, seconds, hours, wall time) and performance.csv are skipped. Numbers are compared exactly.
"""
import os, sys, json, pandas as pd
A, B = sys.argv[1], sys.argv[2]
SKIP_FILES = {'performance.csv', 'orca_runtime.csv'}
TIME_COLS = ('ms', 'seconds', 'hours', 'wall', 'time', 'runtime', 'created')
def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if not any(t in k.lower() for t in TIME_COLS)}
    if isinstance(o, list): return [strip(v) for v in o]
    return o
bad, n = 0, 0
for root, _, files in os.walk(os.path.join(A, 'results')):
    for f in sorted(files):
        if f in SKIP_FILES or f.endswith('.jsonl') and os.path.getsize(os.path.join(root, f)) > 5e8: continue
        pa = os.path.join(root, f); pb = pa.replace(A, B, 1)
        if not os.path.exists(pb): print('MISSING in B:', pb); bad += 1; continue
        n += 1
        if f.endswith('.csv'):
            da, db = pd.read_csv(pa), pd.read_csv(pb)
            cols = [c for c in da.columns if not any(t == c.lower() or c.lower().endswith('_' + t) or c.lower().startswith(t) for t in TIME_COLS)]
            if list(da.columns) != list(db.columns) or not da[cols].equals(db[cols]): print('DIFF', pa); bad += 1
        elif f.endswith('.json'):
            if strip(json.load(open(pa))) != strip(json.load(open(pb))): print('DIFF', pa); bad += 1
        else:
            if open(pa, 'rb').read() != open(pb, 'rb').read(): print('DIFF', pa); bad += 1
print('%d files compared, %d differences' % (n, bad))
