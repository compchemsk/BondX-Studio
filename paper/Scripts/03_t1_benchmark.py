#!/usr/bin/env python3
"""03_t1_benchmark.py
Tier T1 (geometry) of BondX on the HAL59 dimers and on the rigid scans.
For every geometry: BondX geometric score G (bondx_cli.js gscore, same function as the app),
pass/fail of the default criteria, and the materials-path detection (bondx_cli.js detect-xyz).
theta2 = X...Y-R with R the atom bonded to Y that is closest to Y (lone-pair acceptors only).
Output: results/t1_hal59.csv, results/t1_scans.csv
"""
import os, csv, sys
import numpy as np
from bxutil import read_xyz, cli, angle_deg, fragments

def theta2(el, X, xi, yi):
    lab, adj = fragments(el, X, 0.35)
    nb = [j for j in adj[yi] if lab[j] == lab[yi]]
    if not nb: return None
    r = min(nb, key=lambda j: np.linalg.norm(X[j] - X[yi]))
    return angle_deg(X[xi] - X[yi], X[r] - X[yi])

def score_row(xyz, r):
    el, X = read_xyz(xyz); xi = int(r['X']) - 1; ai = int(r['A']) - 1
    if r['ring']:
        ring = [int(i) - 1 for i in str(r['ring']).split(',')]; c = X[ring].mean(axis=0)
        n = np.cross(X[ring[1]] - X[ring[0]], X[ring[2]] - X[ring[0]]); v = X[xi] - c
        alpha = float(np.degrees(np.arccos(abs(v @ n) / np.linalg.norm(v) / np.linalg.norm(n))))
        d = float(np.linalg.norm(v)); t1 = angle_deg(X[ai] - X[xi], c - X[xi])
        g = cli('gscore', '--el', r['halogen'], '--acc', 'AR', '--d', d, '--t1', t1, '--alpha', alpha); t2 = None
    else:
        yi = int(float(r['Y'])) - 1; d = float(np.linalg.norm(X[yi] - X[xi])); t1 = angle_deg(X[ai] - X[xi], X[yi] - X[xi]); t2 = theta2(el, X, xi, yi); alpha = None
        args = ['gscore', '--el', r['halogen'], '--acc', r['acc_el'], '--d', d, '--t1', t1] + (['--t2', t2] if t2 is not None else [])
        g = cli(*args)
    det = cli('detect-xyz', xyz)['hits']
    mat = [h for h in det if h['X'] == int(r['X']) and h['status'] == 'XB']
    return dict(d=round(d, 4), theta1=round(t1, 2), theta2=round(t2, 1) if t2 is not None else '', alpha=round(alpha, 1) if alpha is not None else '',
                covered=g['covered'], cutoff=g['cutoff'], R=round(g['R'], 4) if g['R'] else '', gd=round(g['gd'], 4) if g['gd'] is not None else '',
                gt=round(g['gt'], 4) if g['gt'] is not None else '', g2=round(g['g2'], 4) if g['g2'] is not None else '', pass_T1=g['pass'], G=round(g['G'], 4),
                materials_path_XB=bool(mat))

def run(listfile, folder, outfile, keyname):
    rows = []
    for r in csv.DictReader(open(listfile)):
        fn = os.path.join(folder, r[keyname] if keyname == 'file' else r['name'] + '.xyz')
        res = score_row(fn, r)
        base = {k: r[k] for k in r if k not in ('d', 'theta1', 'alpha')}
        rows.append({**base, **res})
    with open(outfile, 'w', newline='') as f:
        w = csv.DictWriter(f, list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(outfile, len(rows), 'rows')
    return rows

if __name__ == '__main__':
    root = sys.argv[1] if len(sys.argv) > 1 else '..'
    os.makedirs(os.path.join(root, 'results'), exist_ok=True)
    run(os.path.join(root, 'data/hal59/hal59.csv'), os.path.join(root, 'data/hal59'), os.path.join(root, 'results/t1_hal59.csv'), 'name')
    run(os.path.join(root, 'data/scans/scans.csv'), os.path.join(root, 'data/scans'), os.path.join(root, 'results/t1_scans.csv'), 'file')
