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

def guess_hybridization(el, center, nbr_pos, nbr_els):
    n = len(nbr_pos)
    if n >= 2:
        angles = [angle_deg(nbr_pos[i] - center, nbr_pos[j] - center)
                  for i in range(n) for j in range(i + 1, n)]
        avg = sum(angles) / len(angles)
        if n == 2 and avg > 165: return 'sp'
        if n == 3: return 'sp2' if sum(angles) >= 350 else 'sp3'
        if n == 2:
            if el in ('C', 'N'):
                dmin = min(np.linalg.norm(p - center) for p in nbr_pos)
                return 'sp2' if (avg > 114 or dmin < 1.40) else 'sp3'
            return 'sp3'
        return 'sp3'
    if n == 1:
        d = float(np.linalg.norm(nbr_pos[0] - center)); o = nbr_els[0]
        if el == 'O': return 'sp2' if (o == 'C' and d < 1.30) or (o == 'S' and d < 1.50) or (o == 'P' and d < 1.55) or (o == 'N' and d < 1.28) else 'sp3'
        if el == 'N': return 'sp' if d < 1.20 else ('sp2' if d < 1.39 else 'sp3')
        if el == 'C': return 'sp' if d < 1.22 else ('sp2' if d < 1.42 else 'sp3')
        if el == 'S': return 'sp2' if (o == 'C' and d < 1.70) or o == 'P' else 'sp3'
    return 'sp3'

def theta2_and_hyb(el, X, xi, yi):
    lab, adj = fragments(el, X, 0.35)
    nb = [j for j in adj[yi] if lab[j] == lab[yi]]
    if not nb: return None, 'sp2'
    r = min(nb, key=lambda j: np.linalg.norm(X[j] - X[yi]))
    t2 = angle_deg(X[xi] - X[yi], X[r] - X[yi])
    hyb = guess_hybridization(el[yi], X[yi], [X[j] for j in nb], [el[j] for j in nb])
    return t2, hyb

def score_row(xyz, r):
    el, X = read_xyz(xyz); xi = int(r['X']) - 1; ai = int(r['A']) - 1
    donor_is_cbound = el[ai] == 'C'
    if r['ring']:
        ring = [int(i) - 1 for i in str(r['ring']).split(',')]; c = X[ring].mean(axis=0)
        n = np.cross(X[ring[1]] - X[ring[0]], X[ring[2]] - X[ring[0]]); v = X[xi] - c
        alpha = float(np.degrees(np.arccos(abs(v @ n) / np.linalg.norm(v) / np.linalg.norm(n))))
        d = float(np.linalg.norm(v)); t1 = angle_deg(X[ai] - X[xi], c - X[xi])
        g = cli('gscore', '--el', r['halogen'], '--acc', 'AR', '--d', d, '--t1', t1, '--alpha', alpha); t2 = None
    else:
        yi = int(float(r['Y'])) - 1; d = float(np.linalg.norm(X[yi] - X[xi])); t1 = angle_deg(X[ai] - X[xi], X[yi] - X[xi]); t2, hyb = theta2_and_hyb(el, X, xi, yi); alpha = None
        args = ['gscore', '--el', r['halogen'], '--acc', r['acc_el'], '--d', d, '--t1', t1, '--hyb', hyb] + (['--t2', t2] if t2 is not None else [])
        g = cli(*args)
    det = cli('detect-xyz', xyz)['hits']
    mat = [h for h in det if h['X'] == int(r['X']) and h['status'] == 'XB']
    covered = bool(g['covered'] and donor_is_cbound)
    return dict(d=round(d, 4), theta1=round(t1, 2), theta2=round(t2, 1) if t2 is not None else '', theta2_hybridization=hyb if not r['ring'] else '',
                theta2_opt_deg=g['theta2Opt'] if not r['ring'] else '', alpha=round(alpha, 1) if alpha is not None else '', donor_scope_c_bound=donor_is_cbound,
                covered=covered, cutoff=g['cutoff'], R=round(g['R'], 4) if g['R'] else '', gd=round(g['gd'], 4) if g['gd'] is not None else '',
                gt=round(g['gt'], 4) if g['gt'] is not None else '', g2=round(g['g2'], 4) if g['g2'] is not None else '', pass_T1=bool(g['pass'] and donor_is_cbound), G=round(g['G'], 4) if donor_is_cbound else 0,
                materials_path_XB=bool(mat and donor_is_cbound))

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
