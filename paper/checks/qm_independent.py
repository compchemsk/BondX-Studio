#!/usr/bin/env python3
"""Independent check (no BondX code) of the NBO gate sums, EDA terms and Gaussian CP energies in results/*.csv."""
import os, re, json, csv
import pandas as pd
import sys
ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
P = ROOT + '/qm_package'; R = ROOT + '/results'
nbo = pd.read_csv(f'{R}/qm_nbo.csv'); eda = pd.read_csv(f'{R}/qm_eda.csv'); cp = pd.read_csv(f'{R}/qm_cp.csv')
jobs = {(j['set'], j['name'], j['program']): j for f in ('jobs_small.csv', 'jobs_protein.csv') for j in csv.DictReader(open(f'{P}/{f}'))}
SIDE = r'(?:\d+|\*+)\.\s*(BD\*?|LP\*?|CR|RY\*?|LV|3C\*?)\s*\(\s*\d+\s*\)\s*([A-Z][a-z]?)\s*(\d+)(?:\s*-\s*([A-Z][a-z]?)\s*(\d+))?(?:\s*-\s*([A-Z][a-z]?)\s*(\d+))?'
LINE = re.compile(r'^\s*' + SIDE + r'\s*/\s*' + SIDE + r'\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
def entries(txt):
    L = txt.split('\n'); st = max(i for i, l in enumerate(L) if 'SECOND ORDER PERTURBATION THEORY ANALYSIS' in l.upper()); out = []
    for l in L[st:]:
        if 'NATURAL BOND ORBITALS (Summary)' in l: break
        m = LINE.match(l)
        if m:
            g = m.groups()
            dat = (g[0], [int(x) for x in (g[2], g[4], g[6]) if x]); acc = (g[7], [int(x) for x in (g[9], g[11], g[13]) if x])
            out.append((dat, acc, float(g[14])))
    return out
bad = 0; n = 0
for (s, nm), grp in nbo.groupby(['set', 'name']):
    j = jobs[(s, nm, 'gaussian')]; m = json.load(open(f"{P}/{j['manifest']}")); txt = open(f"{P}/{os.path.splitext(j['input'])[0]}.log", errors='replace').read(); E = entries(txt)
    for xb in m['xbs']:
        X, A = xb['X'], xb['A']; ring = set(xb.get('ring') or [])
        has_sig = any(a[0] == 'BD*' and sorted(a[1]) == sorted([X, A]) for _, a, _ in E)
        tot = 0.0
        for d, a, e in E:
            to_sig = (a[0] == 'BD*' and sorted(a[1]) == sorted([X, A])) or (not has_sig and a[0] == 'LP*' and a[1] == [X])
            to_ry = a[0].startswith('RY') and a[1] == [X]
            if not (to_sig or to_ry): continue
            if xb['type'] == 'atom' and d[0] == 'LP' and d[1][0] == xb['Y']: tot += e
            elif xb['type'] == 'pi' and d[0] == 'BD' and len(d[1]) == 2 and set(d[1]) <= ring and to_sig: tot += e
        row = grp[grp.xb == xb['id']]; n += 1
        if len(row) != 1 or abs(row.E2_core.iloc[0] - round(tot, 3)) > 0.006: bad += 1; print('NBO mismatch', s, nm, xb['id'], tot, row.E2_core.tolist())
print('NBO gate sums checked:', n, 'mismatches:', bad)
b2 = 0
for r in eda.itertuples():
    j = jobs[(r.set, r.name, 'qchem_eda')]; t = open(f"{P}/{os.path.splitext(j['input'])[0]}.out", errors='replace').read()
    i = t.rfind('Simplified EDA Summary'); blk = t[i:i + 1500]
    v = {k: float(re.search(k + r'\s+(-?\d+\.\d+)', blk).group(1)) for k in ('FROZEN', 'POLARIZATION', 'CHARGE TRANSFER', 'TOTAL')}
    if abs(v['FROZEN'] - r.frz) > 1e-3 or abs(v['POLARIZATION'] - r.pol) > 1e-3 or abs(v['CHARGE TRANSFER'] - r.ct) > 1e-3 or abs(v['TOTAL'] - r.total) > 1e-3: b2 += 1; print('EDA mismatch', r.name)
    if abs(v['FROZEN'] + v['POLARIZATION'] + v['CHARGE TRANSFER'] - v['TOTAL']) > 2e-3: print('EDA sum?', r.name)
print('EDA rows checked:', len(eda), 'mismatches:', b2)
b3 = 0; g = cp[cp.program == 'gaussian']
for r in g.itertuples():
    j = jobs[(r.set, r.name, 'gaussian_cp')]; t = open(f"{P}/{os.path.splitext(j['input'])[0]}.log", errors='replace').read()
    k = float(re.findall(r'complexation energy =\s+(-?\d+\.\d+) kcal/mole \(corrected\)', t)[-1])
    if abs(k - r.Eint_CP_kJ / 4.184) > 0.006: b3 += 1; print('CP mismatch', r.name, k, r.Eint_CP_kJ / 4.184)
print('Gaussian CP rows checked:', len(g), 'mismatches:', b3)
