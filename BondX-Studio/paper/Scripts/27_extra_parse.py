#!/usr/bin/env python3
"""27_extra_parse.py
Independent parser (no BondX code) for the additional Gaussian 16 NBO jobs (sets G, J, and K).
For every output it reads the geometry and the second-order perturbation table of NBO 3.1 and computes
  E2_valence : sum of LP(Y) -> BD*(A-X) or LP*(X) terms (Y = T1 acceptor atom, X = halogen, A = its carbon)
  E2_RY      : sum of LP(Y) -> RY*(X)
  hb_Y       : sum of LP(X) -> BD*(Y-H) terms (competing hydrogen bond from the acceptor atom)
  hb_any     : sum of LP(X) -> BD*(N-H or O-H) terms to atoms outside the halogen-bearing molecule
  LP_on_Y    : NBO prints at least one lone pair on Y
Sets G and J take X, A, and Y from the BondX manifests. Set K (published hIYD inputs with a new acetyl cap) has no
manifest: X is the iodine, A its carbon, and Y the backbone N nearest to X.
Usage: python3 27_extra_parse.py ../qm_extra ../results/extra
Output: results/extra/extra_parsed.csv
"""
import os, re, sys, json, glob, math
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else '../qm_extra'
OUT = sys.argv[2] if len(sys.argv) > 2 else '../results/extra'
EL = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I'}
ROW = re.compile(r'^\s*(?:\d+|\*+)\.\s+(\S+)\s*\(\s*\d+\)\s*(.+?)\s+/\s*(?:\d+|\*+)\.\s+(\S+)\s*\(\s*\d+\)\s*(.+?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$')
ATOMS = re.compile(r'([A-Z][a-z]?)\s*(\d+)')


def geometry(lines):
    """Last 'Input orientation' block (nosymm jobs)."""
    start = None
    for i, l in enumerate(lines):
        if 'Input orientation:' in l or 'Standard orientation:' in l:
            start = i
    xyz = []
    for l in lines[start + 5:]:
        if l.strip().startswith('----'):
            break
        p = l.split()
        xyz.append((EL.get(int(p[1]), p[1]), float(p[3]), float(p[4]), float(p[5])))
    return xyz


def table(lines):
    """Second-order table rows: (donor type, donor atoms, acceptor type, acceptor atoms, E2)."""
    rows, on = [], False
    for l in lines:
        if 'Second Order Perturbation Theory Analysis' in l:
            on = True; rows = []; continue
        if on and 'Natural Bond Orbitals (Summary)' in l:
            on = False; continue
        if on:
            m = ROW.match(l)
            if m:
                da = [int(n) for _, n in ATOMS.findall(m.group(2))]
                aa = [int(n) for _, n in ATOMS.findall(m.group(4))]
                rows.append((m.group(1), da, m.group(3), aa, float(m.group(5))))
    return rows


def lone_pairs(lines):
    """Atoms that carry a Lewis lone pair in the NBO listing, e.g. '115. (1.99308) LP ( 1) I  23'."""
    s = set()
    for l in lines:
        m = re.match(r'^\s*(?:\d+|\*+)\.\s+\(\s*[0-9.]+\)\s+LP\s*\(\s*\d+\)\s*[A-Z][a-z]?\s*(\d+)', l)
        if m:
            s.add(int(m.group(1)))
    return s


def dist(a, b):
    return math.dist(a[1:], b[1:])


def bonded(xyz, i, j):
    lim = {'H': 0.32, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'P': 1.07, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39, 'F': 0.57}
    return dist(xyz[i - 1], xyz[j - 1]) < lim.get(xyz[i - 1][0], 0.8) + lim.get(xyz[j - 1][0], 0.8) + 0.45


def molecule(xyz, start):
    seen, stack = {start}, [start]
    while stack:
        a = stack.pop()
        for b in range(1, len(xyz) + 1):
            if b not in seen and bonded(xyz, a, b):
                seen.add(b); stack.append(b)
    return seen


def analyse(log, X, A, Y):
    L = open(log, errors='ignore').read().split('\n')
    xyz = geometry(L); T = table(L); lp = lone_pairs(L)
    mol1 = molecule(xyz, X)
    val = ry = hbY = hbany = 0.0
    for dt, da, at, aa, e in T:
        if dt == 'LP' and da == [Y]:
            if at == 'BD*' and sorted(aa) == sorted([A, X]):
                val += e
            elif at == 'LP*' and aa == [X]:
                val += e
            elif at == 'RY*' and aa == [X]:
                ry += e
        if dt == 'LP' and da == [X] and at == 'BD*' and len(aa) == 2:
            heavy = [a for a in aa if xyz[a - 1][0] != 'H']; hyd = [a for a in aa if xyz[a - 1][0] == 'H']
            if len(heavy) == 1 and len(hyd) == 1 and xyz[heavy[0] - 1][0] in ('N', 'O') and heavy[0] not in mol1:
                hbany += e
                if heavy[0] == Y:
                    hbY += e
    ok = any('Normal termination' in l for l in L[-5:])
    d = dist(xyz[X - 1], xyz[Y - 1])
    v1 = [xyz[A - 1][k] - xyz[X - 1][k] for k in (1, 2, 3)]; v2 = [xyz[Y - 1][k] - xyz[X - 1][k] for k in (1, 2, 3)]
    th = math.degrees(math.acos(sum(p * q for p, q in zip(v1, v2)) / (math.hypot(*v1) * math.hypot(*v2))))
    nheavy = sum(1 for b in range(1, len(xyz) + 1) if b != Y and bonded(xyz, Y, b) and xyz[b - 1][0] != 'H')
    nh = sum(1 for b in range(1, len(xyz) + 1) if b != Y and bonded(xyz, Y, b) and xyz[b - 1][0] == 'H')
    return dict(normal=ok, nNBOrows=len(T), E2_valence=round(val, 2), E2_RY=round(ry, 2), hb_Y=round(hbY, 2), hb_any=round(hbany, 2),
                LP_on_Y=Y in lp, d=round(d, 3), theta1=round(th, 1), Y_el=xyz[Y - 1][0], Y_heavy_nb=nheavy, Y_H=nh, natoms=len(xyz))


rows = []
for s in ('G_hard_negatives', 'J_4dst_amide'):
    for man in sorted(glob.glob(os.path.join(ROOT, s, 'manifests', '*_manifest.json'))):
        m = json.load(open(man)); name = os.path.basename(man).replace('_manifest.json', '')
        log = os.path.join(ROOT, s, 'gaussian', name + '.log')
        if not os.path.exists(log):
            continue
        xb = m['xbs'][0]
        r = analyse(log, xb['X'], xb['A'], xb['Y']); r.update(set=s, name=name, contact=xb['label'], charge=m['totalCharge'])
        rows.append(r)
for log in sorted(glob.glob(os.path.join(ROOT, 'K_hiyd_amide', 'gaussian', '*.log'))):
    L = open(log, errors='ignore').read().split('\n'); xyz = geometry(L)
    X = next(i + 1 for i, a in enumerate(xyz) if a[0] == 'I')
    A = min((i + 1 for i, a in enumerate(xyz) if a[0] == 'C'), key=lambda j: dist(xyz[X - 1], xyz[j - 1]))
    cand = [i + 1 for i, a in enumerate(xyz) if a[0] == 'N' and (i + 1) not in molecule(xyz, X)]
    Y = min(cand, key=lambda j: dist(xyz[X - 1], xyz[j - 1]))
    r = analyse(log, X, A, Y); r.update(set='K_hiyd_amide', name=os.path.basename(log)[:-4], contact='I···N (nearest backbone N)', charge=None)
    rows.append(r)
df = pd.DataFrame(rows)
os.makedirs(OUT, exist_ok=True)
df.to_csv(os.path.join(OUT, 'extra_parsed.csv'), index=False)
print(df[['set', 'name', 'normal', 'E2_valence', 'E2_RY', 'hb_Y', 'hb_any', 'LP_on_Y', 'd', 'theta1', 'Y_el', 'Y_heavy_nb', 'Y_H']].to_string())
