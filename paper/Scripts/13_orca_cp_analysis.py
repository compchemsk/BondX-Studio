#!/usr/bin/env python3
"""13_orca_cp_analysis.py
Analysis of the ORCA 5.0.1 counterpoise jobs (wB97X-V/def2-TZVPD, RIJCOSX def2/J, gas phase).
E_int = E_AB - E_A(AB basis) - E_B(AB basis). Binding energy BE = -E_int (positive = bound).

Reads   results/qm_cp.csv (07_collect_results.py), data/hal59/hal59.csv, results/t1_hal59.csv,
        results/t1_scans.csv, results/protein_models.csv, results/qm_status.csv and the ORCA outputs.
Writes  results/orca_cp_hal59.csv, results/orca_cp_scans.csv, results/orca_cp_protein.csv,
        results/orca_runtime.csv, results/orca_cp_summary.json
Usage   python3 13_orca_cp_analysis.py ..
"""
import os, sys, csv, json, re, glob
import numpy as np
from scipy import stats

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
RES = os.path.join(ROOT, 'results'); PKG = os.path.join(ROOT, 'qm_package')
KJ_PER_KCAL = 4.184


def rd(p): return list(csv.DictReader(open(p)))


def wr(p, rows):
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)


def errstats(x, y):
    """x = method, y = reference (kcal/mol)."""
    x, y = np.asarray(x, float), np.asarray(y, float); e = x - y
    return dict(n=int(len(e)), MSE=float(e.mean()), MAE=float(np.abs(e).mean()), RMSD=float(np.sqrt((e ** 2).mean())),
                MaxAE=float(np.abs(e).max()), pearson_r=float(stats.pearsonr(x, y)[0]) if len(e) > 2 else None,
                spearman_rho=float(stats.spearmanr(x, y)[0]) if len(e) > 2 else None,
                MRAE_pct=float(np.mean(np.abs(e) / np.abs(y)) * 100))


cp = {(r['set'], r['name']): float(r['Eint_CP_kJ']) / KJ_PER_KCAL for r in rd(os.path.join(RES, 'qm_cp.csv')) if r['program'] == 'orca'}
summary = {'level': 'wB97X-V/def2-TZVPD (ORCA 5.0.1, RIJCOSX def2/J, TightSCF, gas phase), counterpoise corrected', 'n_systems': len(cp)}

# ---------------- HAL59 ----------------
hal = {r['name']: r for r in rd(os.path.join(ROOT, 'data', 'hal59', 'hal59.csv'))}
t1 = {r['name']: r for r in rd(os.path.join(RES, 't1_hal59.csv'))}
H = []
for n, r in hal.items():
    e = cp[('A_HAL59', n)]
    H.append(dict(name=n, subset=r['subset'], role=r['role'], halogen=r['halogen'].capitalize(), acc=r['acc_el'] if r['acc_el'] else 'ring',
                  d=r['d'], theta1=r['theta1'], ref_kcal=float(r['ref_kcal']), Eint_kcal=round(e, 3), BE_kcal=round(-e, 3),
                  err_kcal=round(-e - float(r['ref_kcal']), 3), G=t1[n]['G'], covered=t1[n]['covered'], pass_T1=t1[n]['pass_T1']))
wr(os.path.join(RES, 'orca_cp_hal59.csv'), H)
be = [h['BE_kcal'] for h in H]; ref = [h['ref_kcal'] for h in H]
summary['hal59_all'] = errstats(be, ref)
for key, sel in (('XB51', lambda h: h['subset'] == 'XB51'), ('X40', lambda h: h['subset'] != 'XB51'),
                 ('Cl', lambda h: h['halogen'] == 'Cl'), ('Br', lambda h: h['halogen'] == 'Br'), ('I', lambda h: h['halogen'] == 'I'),
                 ('XB_dimers', lambda h: h['role'] == 'halogen bond'), ('controls', lambda h: h['role'] != 'halogen bond'),
                 ('covered', lambda h: h['covered'] == 'True')):
    s = [h for h in H if sel(h)]
    if len(s) >= 2: summary['hal59_' + key] = errstats([h['BE_kcal'] for h in s], [h['ref_kcal'] for h in s])
summary['hal59_subsets'] = {k: sum(1 for h in H if h['subset'] == k) for k in sorted({h['subset'] for h in H})}
summary['hal59_roles'] = {k: sum(1 for h in H if h['role'] == k) for k in sorted({h['role'] for h in H})}
big = sorted(H, key=lambda h: -abs(h['err_kcal']))[:5]
summary['hal59_largest_errors'] = [(h['name'], h['ref_kcal'], h['BE_kcal'], h['err_kcal']) for h in big]
cov = [h for h in H if h['covered'] == 'True']
summary['G_vs_BE_spearman_covered'] = float(stats.spearmanr([float(h['G']) for h in cov], [h['BE_kcal'] for h in cov])[0])
summary['G_vs_ref_spearman_covered'] = float(stats.spearmanr([float(h['G']) for h in cov], [h['ref_kcal'] for h in cov])[0])
summary['n_covered'] = len(cov)
summary['controls'] = [(h['name'], h['ref_kcal'], h['BE_kcal'], h['pass_T1']) for h in H if h['role'] != 'halogen bond']

# ---------------- scans ----------------
sc = rd(os.path.join(RES, 't1_scans.csv')); S = []
for r in sc:
    st = 'B_scan_distance' if r['scan'] == 'distance' else 'C_scan_angle'
    nm = r['file'][:-4]; e = cp[(st, nm)]
    S.append(dict(base=r['base'], scan=r['scan'], s=r['s'], theta=r['theta_target'], d=r['d'], theta1=r['theta1'], halogen=r['halogen'].capitalize(),
                  acc=r['acc_el'] or 'ring', Eint_kcal=round(e, 3), BE_kcal=round(-e, 3), G=r['G'], pass_T1=r['pass_T1']))
# relative values
for base in {x['base'] for x in S}:
    for scan, key, ref0 in (('distance', 's', '1.0'), ('angle', 'theta', '180')):
        grp = [x for x in S if x['base'] == base and x['scan'] == scan]
        if not grp: continue
        r0 = [x for x in grp if float(x[key]) == float(ref0)][0]['BE_kcal']
        for x in grp: x['BE_rel'] = round(x['BE_kcal'] / r0, 4) if r0 else None
wr(os.path.join(RES, 'orca_cp_scans.csv'), S)
# consistency: s = 1.00 equals the HAL59 geometry
dd = [abs(x['BE_kcal'] - [h['BE_kcal'] for h in H if h['name'] == x['base']][0]) for x in S if x['scan'] == 'distance' and float(x['s']) == 1.0]
summary['scan_s100_vs_hal59_maxdiff_kcal'] = max(dd)
dist = {}
for base in sorted({x['base'] for x in S if x['scan'] == 'distance'}):
    g = sorted([x for x in S if x['base'] == base and x['scan'] == 'distance'], key=lambda x: float(x['s']))
    smin = max(g, key=lambda x: x['BE_kcal'])
    dist[base] = dict(s_at_max_BE=float(smin['s']), BE_max=smin['BE_kcal'], BE_s100=[x['BE_kcal'] for x in g if float(x['s']) == 1.0][0],
                      BE_s090=[x['BE_kcal'] for x in g if float(x['s']) == 0.9][0], BE_s150=[x['BE_kcal'] for x in g if float(x['s']) == 1.5][0],
                      frac_s120=[x['BE_rel'] for x in g if float(x['s']) == 1.2][0], frac_s150=[x['BE_rel'] for x in g if float(x['s']) == 1.5][0],
                      G_zero_from_s=min([float(x['s']) for x in g if float(x['G']) == 0] or [None]))
summary['distance_scans'] = dist
ang = {}
for base in sorted({x['base'] for x in S if x['scan'] == 'angle'}):
    g = sorted([x for x in S if x['base'] == base and x['scan'] == 'angle'], key=lambda x: -float(x['theta']))
    ang[base] = dict(BE_180=g[0]['BE_kcal'], **{f'frac_{int(float(x["theta"]))}': x['BE_rel'] for x in g},
                     monotonic=all(g[i]['BE_kcal'] >= g[i + 1]['BE_kcal'] - 1e-6 for i in range(len(g) - 1)))
summary['angle_scans'] = ang
f140 = [v['frac_140'] for v in ang.values()]; f120 = [v['frac_120'] for v in ang.values()]
summary['angle_frac140_range'] = [min(f140), max(f140)]; summary['angle_frac120_range'] = [min(f120), max(f120)]
# Spearman between G and BE along each distance scan, for points inside the window (G > 0)
rhos = []
for base in dist:
    g = [x for x in S if x['base'] == base and x['scan'] == 'distance']
    rhos.append(float(stats.spearmanr([float(x['G']) for x in g], [x['BE_kcal'] for x in g])[0]))
summary['distance_scan_G_BE_spearman_range'] = [min(rhos), max(rhos)]

# ---------------- protein clusters ----------------
pm = {r['model']: r for r in rd(os.path.join(RES, 'protein_models.csv'))}
P = []
for (st, nm), e in sorted(cp.items()):
    if st != 'F_protein': continue
    m = pm[nm]
    P.append(dict(model=nm, system=m['system'], ligand=m['ligand'], n_atoms=int(m['n_atoms']), charge_donor=m['charge_donor'], charge_acceptor=m['charge_acceptor'],
                  n_xb=m['n_xb'], xbs=m['xbs'], Eint_kJ=round(e * KJ_PER_KCAL, 2), Eint_kcal=round(e, 2)))
if P: wr(os.path.join(RES, 'orca_cp_protein.csv'), P)
elif os.path.exists(os.path.join(RES, 'orca_cp_protein.csv')): os.remove(os.path.join(RES, 'orca_cp_protein.csv'))
summary['protein'] = {p['model']: p['Eint_kcal'] for p in P} if P else 'none'
if P:
    rep = {}
    for a_, b_ in (('1p5e__TBSA__cluster_pcm4', '1p5e__TBSC__cluster_pcm4'), ('2w0s__BVPA__cluster_pcm4', '2w0s__BVPB__cluster_pcm4'), ('4agl__P84A__cluster_pcm4', '4agl__P84B__cluster_pcm4')):
        d = {p['model']: p['Eint_kcal'] for p in P}
        if a_ in d and b_ in d: rep[a_.split('__')[0]] = round(abs(d[a_] - d[b_]), 2)
    summary['protein_replica_absdiff_kcal'] = rep

# ---------------- run times ----------------
T = []
for r in rd(os.path.join(RES, 'qm_status.csv')):
    if r['program'] != 'orca_cp' or r['status'] != 'ok': continue
    out = os.path.join(PKG, os.path.splitext(r['input'])[0] + '.out'); t = open(out, errors='replace').read()
    m = re.search(r'TOTAL RUN TIME:\s*(\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec', t)
    sec = int(m[1]) * 86400 + int(m[2]) * 3600 + int(m[3]) * 60 + int(m[4]) + int(m[5]) / 1000
    nb = re.search(r'Number of basis functions\s+\.\.\.\s+(\d+)', t)
    it = open(os.path.join(PKG, r['input'])).read(); mm = re.search(r'^\* xyz.*$', it, re.M)
    nat = len(it[mm.end():it.index('\n*', mm.end())].strip().split('\n'))
    T.append(dict(set=r['set'], input=r['input'], atoms=nat, basis_functions=int(nb[1]) if nb else '', seconds=round(sec, 1)))
wr(os.path.join(RES, 'orca_runtime.csv'), T)
tot = sum(x['seconds'] for x in T)
summary['runtime'] = dict(jobs=len(T), wall_hours=round(tot / 3600, 1), core_hours_40=round(tot / 3600 * 40, 0),
                          by_set={st: dict(jobs=sum(1 for x in T if x['set'] == st), wall_h=round(sum(x['seconds'] for x in T if x['set'] == st) / 3600, 2),
                                           median_s=float(np.median([x['seconds'] for x in T if x['set'] == st])),
                                           max_s=max(x['seconds'] for x in T if x['set'] == st)) for st in sorted({x['set'] for x in T})},
                          atoms_range=[min(x['atoms'] for x in T), max(x['atoms'] for x in T)],
                          nbf_range=[min(x['basis_functions'] for x in T), max(x['basis_functions'] for x in T)])
json.dump(summary, open(os.path.join(RES, 'orca_cp_summary.json'), 'w'), indent=1)
print(json.dumps(summary, indent=1))
