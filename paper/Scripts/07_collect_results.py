#!/usr/bin/env python3
"""07_collect_results.py
Reads every finished QM output in qm_package/ with the BondX parsers (bondx_cli.js gate / eda) and writes
  results/qm_nbo.csv  one row per NBO job and halogen bond: SumE(2) terms, gate decision, NPA charges, checks
  results/qm_eda.csv  one row per ALMO-EDA job: FRZ, ELEC, PAULI, DISP, POL, CT, TOTAL (kJ/mol), chi_CT, mechanism
  results/qm_cp.csv   counterpoise-corrected interaction energies from Gaussian (Counterpoise=2) and ORCA (AB, A, B)
  results/qm_status.csv  every job with its status (ok / failed / missing)
Outputs are expected next to the inputs: Q-Chem and ORCA *.out, Gaussian *.log.
Usage: python3 07_collect_results.py ../qm_package ../results [--threshold 0.5]
"""
import os, sys, csv, json
from bxutil import cli

PKG = sys.argv[1] if len(sys.argv) > 1 else '../qm_package'
RES = sys.argv[2] if len(sys.argv) > 2 else '../results'
THR = sys.argv[sys.argv.index('--threshold') + 1] if '--threshold' in sys.argv else '0.5'
OK_MARK = {'qchem': 'Thank you very much for using Q-Chem', 'qchem_eda': 'Thank you very much for using Q-Chem',
           'gaussian': 'Normal termination of Gaussian', 'orca': 'ORCA TERMINATED NORMALLY',
           'gaussian_cp': 'Normal termination of Gaussian', 'orca_cp': 'ORCA TERMINATED NORMALLY'}
HARTREE_KJ = 2625.49964


def out_of(inp, prog):
    stem = os.path.splitext(inp)[0]
    return stem + ('.log' if prog.startswith('gaussian') else '.out')


def main():
    jobs = []
    for f in ('jobs_small.csv', 'jobs_protein.csv'):
        p = os.path.join(PKG, f)
        if os.path.exists(p): jobs += list(csv.DictReader(open(p)))
    import re
    nbo_rows, eda_rows, status, cp_rows, orca_e = [], [], [], [], {}
    for j in jobs:
        out = os.path.join(PKG, out_of(j['input'], j['program']))
        if not os.path.exists(out):
            status.append(dict(j, status='missing')); continue
        txt = open(out, errors='replace').read()
        finished = OK_MARK[j['program']] in txt
        if j['program'] == 'gaussian_cp':
            # Gaussian 16 (Counterpoise=2) prints the CP-corrected total energy and the sum of the fragment energies in their own basis
            mc = re.findall(r'Counterpoise corrected energy =\s+(-?\d+\.\d+)', txt); mf = re.findall(r'sum of fragments =\s+(-?\d+\.\d+)', txt)
            mk = re.findall(r'complexation energy =\s+(-?\d+\.\d+) kcal/mole \(corrected\)', txt)
            if not (mc and mf): status.append(dict(j, status='failed: no counterpoise energy')); continue
            eint = (float(mc[-1]) - float(mf[-1])) * HARTREE_KJ
            if mk: assert abs(eint / 4.184 - float(mk[-1])) < 0.006, (j['input'], eint / 4.184, mk[-1])   # agrees with the printed kcal/mol value
            cp_rows.append(dict(set=j['set'], name=j['name'], program='gaussian', level='wB97XD/def2TZVP', Eint_CP_kJ=round(eint, 3), finished=finished))
            status.append(dict(j, status='ok' if finished else 'not finished')); continue
        if j['program'] == 'orca_cp':
            m = re.findall(r'FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)', txt)
            if not m: status.append(dict(j, status='failed: no final energy')); continue
            part = re.search(r'_cp_(AB|A|B)\.inp$', j['input']).group(1); orca_e.setdefault((j['set'], j['name']), {})[part] = float(m[-1])
            status.append(dict(j, status='ok' if finished else 'not finished')); continue
        if j['program'] == 'qchem_eda':
            r = cli('eda', out)
            if r.get('error'):
                status.append(dict(j, status='failed: ' + r['error'][:80])); continue
            d = r.get('descriptors') or {}
            eda_rows.append(dict(set=j['set'], name=j['name'], finished=finished, units=r.get('units'), frz=r.get('frz'), elec=r.get('elec'), pauli=r.get('pauli'), elec_pauli=r.get('elecPauli'), disp=r.get('disp'),
                                 pol=r.get('pol'), ct=r.get('ct'), total=r.get('total'), chi_CT=d.get('chiCT'), ct_pol_ratio=d.get('ctPolRatio'), mechanism=d.get('mechanism')))
        else:
            r = cli('gate', os.path.join(PKG, j['manifest']), out, '--threshold', THR)
            if r.get('error'):
                status.append(dict(j, status='failed: ' + r['error'][:80])); continue
            v = r['verification']
            for g in r['gate']:
                nbo_rows.append(dict(set=j['set'], name=j['name'], program=j['program'], finished=finished, xb=g['id'], label=g['label'], E2_core=round(g['E2_core'], 3),
                                     E2_total=round(g['E2_total'], 3), lp_sigma=round(g['contrib']['lp_sigma'], 3), lp_ry=round(g['contrib']['lp_ry'], 3), pi_sigma=round(g['contrib']['pi_sigma'], 3),
                                     verified=g['verified'], ionicX=g.get('ionicX'), npaX=g['npaX'], npaY=g['npaY'], ct_NPA=(r['fragNPA'] or {}).get('ct'), print_threshold=r['printThreshold'],
                                     geometry_ok=v.get('ok'), geometry_note=v.get('note')))
        status.append(dict(j, status='ok' if finished else 'not finished'))
    for (st, nm), e in orca_e.items():
        if all(k in e for k in ('AB', 'A', 'B')):
            cp_rows.append(dict(set=st, name=nm, program='orca', level='wB97X-V/def2-TZVPD', Eint_CP_kJ=round((e['AB'] - e['A'] - e['B']) * HARTREE_KJ, 3), finished=True))
    os.makedirs(RES, exist_ok=True)
    for fn, rows in (('qm_nbo.csv', nbo_rows), ('qm_eda.csv', eda_rows), ('qm_cp.csv', cp_rows), ('qm_status.csv', status)):
        if not rows: print('no rows for', fn); continue
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(os.path.join(RES, fn), 'w', newline='') as f:
            w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
        print(fn, len(rows))
    from collections import Counter
    print(Counter((s['set'], s['program'], s['status'].split(':')[0]) for s in status))


if __name__ == '__main__':
    main()
