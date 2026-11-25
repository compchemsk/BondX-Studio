#!/usr/bin/env python3
"""02_make_scans.py
Rigid scans of halogen-bond geometry for a subset of HAL59 dimers.
  distance scan: the acceptor molecule moves along the X...Y line; d = s * d_eq, s = 0.90 ... 1.50 (8 points)
  angle scan   : the acceptor molecule rotates about X; theta1(A-X...Y) = 180, 170, ..., 120 deg (7 points).
                 The rotation axis lies in the donor ring plane and is perpendicular to A-X, so the acceptor
                 orientation relative to the X...Y line (theta2) does not change.
Output: data/scans/<name>_d<s>.xyz, data/scans/<name>_a<theta>.xyz and data/scans/scans.csv
Usage : python3 02_make_scans.py data/hal59 data/scans
"""
import os, sys, csv
import numpy as np
from bxutil import read_xyz, analyse, write_xyz, fragments, angle_deg

DIST_SET = ['19_benCl-acetone', '20_benBr-acetone', '21_benI-acetone', '22_benCl-Nm3', '23_benBr-Nm3', '24_benI-Nm3',
            '26_benI-mSH', '28_CH3I-ben', 'NH3_PhI', 'MeI_pyr', 'FI_pyr', 'NH3_FBr']
ANGLE_SET = ['19_benCl-acetone', '20_benBr-acetone', '21_benI-acetone', '22_benCl-Nm3', '23_benBr-Nm3', '24_benI-Nm3']
SCALES = [0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.35, 1.50]
ANGLES = [180, 170, 160, 150, 140, 130, 120]

def rot(axis, phi):
    a = axis / np.linalg.norm(axis); K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(phi) * K + (1 - np.cos(phi)) * K @ K

def split(el, X, info):
    lab = None
    for tol in (0.45, 0.35, 0.25, 0.15, 0.05):
        lab, _ = fragments(el, X, tol)
        if len(set(lab)) == 2 and lab.count(lab[int(info['X']) - 1]) == int(info['frag_donor']): break
    don = [i for i in range(len(el)) if lab[i] == lab[int(info['X']) - 1]]
    acc = [i for i in range(len(el)) if lab[i] != lab[int(info['X']) - 1]]
    return don, acc

def target(X, info):
    if info['ring']: return X[[int(i) - 1 for i in str(info['ring']).split(',')]].mean(axis=0)
    return X[int(info['Y']) - 1]

def min_contact(X, don, acc):
    return float(min(np.linalg.norm(X[i] - X[j]) for i in don for j in acc))

def main(src, dst):
    os.makedirs(dst, exist_ok=True); rows = []
    meta = {r['name']: r for r in csv.DictReader(open(os.path.join(src, 'hal59.csv')))}
    for name in DIST_SET:
        el, X0 = read_xyz(os.path.join(src, name + '.xyz')); info = meta[name]
        don, acc = split(el, X0, info); x = X0[int(info['X']) - 1]; y = target(X0, info)
        u = (y - x) / np.linalg.norm(y - x); d0 = np.linalg.norm(y - x)
        for s in SCALES:
            X = X0.copy(); X[acc] += u * (s - 1.0) * d0
            fn = '%s_d%03d.xyz' % (name, round(s * 100)); write_xyz(os.path.join(dst, fn), el, X, '%s distance scan s=%.2f' % (name, s))
            a = dict(d=float(np.linalg.norm(target(X, info) - X[int(info['X']) - 1])), theta1=angle_deg(X[int(info['A']) - 1] - X[int(info['X']) - 1], target(X, info) - X[int(info['X']) - 1]))
            rows.append(dict(file=fn, base=name, scan='distance', s=s, theta_target='', d=round(a['d'], 4), theta1=round(a['theta1'], 2),
                             min_contact=round(min_contact(X, don, acc), 3), X=info['X'], A=info['A'], Y=info['Y'], ring=info['ring'], halogen=info['halogen'], acc_el=info['acc_el']))
    for name in ANGLE_SET:
        el, X0 = read_xyz(os.path.join(src, name + '.xyz')); info = meta[name]
        don, acc = split(el, X0, info); xi, ai = int(info['X']) - 1, int(info['A']) - 1; x = X0[xi]
        ring = [i for i in don if el[i] == 'C'][:6]
        nrm = np.cross(X0[ring[1]] - X0[ring[0]], X0[ring[2]] - X0[ring[0]]); ax = np.cross(X0[ai] - x, nrm)
        th_eq = angle_deg(X0[ai] - x, target(X0, info) - x)
        for th in ANGLES:
            R = rot(ax, np.radians(th_eq - th)); X = X0.copy(); X[acc] = (X0[acc] - x) @ R.T + x
            th_now = angle_deg(X[ai] - x, target(X, info) - x)
            if abs(th_now - th) > 0.5:      # rotate the other way if needed
                R = rot(ax, -np.radians(th_eq - th)); X = X0.copy(); X[acc] = (X0[acc] - x) @ R.T + x; th_now = angle_deg(X[ai] - x, target(X, info) - x)
            fn = '%s_a%03d.xyz' % (name, th); write_xyz(os.path.join(dst, fn), el, X, '%s angle scan theta1=%d' % (name, th))
            rows.append(dict(file=fn, base=name, scan='angle', s=1.0, theta_target=th, d=round(float(np.linalg.norm(target(X, info) - x)), 4), theta1=round(th_now, 2),
                             min_contact=round(min_contact(X, don, acc), 3), X=info['X'], A=info['A'], Y=info['Y'], ring=info['ring'], halogen=info['halogen'], acc_el=info['acc_el']))
    with open(os.path.join(dst, 'scans.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print('%d scan geometries written to %s' % (len(rows), dst))

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'data/hal59', sys.argv[2] if len(sys.argv) > 2 else 'data/scans')
