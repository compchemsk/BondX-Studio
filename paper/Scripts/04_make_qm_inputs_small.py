#!/usr/bin/env python3
"""04_make_qm_inputs_small.py
QM inputs for the gas-phase benchmark. BondX builds every model (bondx_cli.js build-xyz).
  A_HAL59          59 dimers   NBO  wB97X-D/def2-TZVP   Gaussian 16, ORCA, Q-Chem
                               EDA  wB97X-V/def2-TZVPD  Q-Chem (ALMO-EDA2 with EDA_BSSE)
                               with every EDA job: counterpoise interaction energy in Gaussian 16 (wB97XD/def2TZVP)
                               and ORCA (wB97X-V/def2-TZVPD; three single points AB, A(AB basis), B(AB basis))
  B_scan_distance  96 points   NBO (3 programs) + EDA (Q-Chem), same levels
  C_scan_angle     42 points   NBO (3 programs) + EDA (Q-Chem), same levels
  D_method_basis   12 dimers x {wB97X-D, M06-2X, B3LYP-D3BJ, PBE0-D3BJ, HF}/def2-TZVP
                               + wB97X-D/{def2-SVP, def2-TZVPD, def2-QZVP}           NBO, 3 programs
                               (Gaussian 16 has no built-in def2-TZVPD: ORCA and Q-Chem only)
  E_eda_settings   12 dimers   EDA wB97X-V/{def2-SVPD, def2-QZVP} and def2-TZVPD without BSSE (Q-Chem)
Usage: python3 04_make_qm_inputs_small.py ../data ../qm_package
"""
import os, sys, csv
from collections import Counter
from bxutil import cli

SENS = ['19_benCl-acetone', '20_benBr-acetone', '21_benI-acetone', '22_benCl-Nm3', '23_benBr-Nm3', '24_benI-Nm3',
        '26_benI-mSH', '28_CH3I-ben', 'NH3_PhI', 'MeI_pyr', 'FI_pyr', 'NH3_FBr']
METHODS = ['wB97X-D', 'M06-2X', 'B3LYP-D3BJ', 'PBE0-D3BJ', 'HF']
BASES = ['def2-SVP', 'def2-TZVPD', 'def2-QZVP']
RES = dict(nproc=40, mem=50000)
jobs = []


def atomargs(r):
    a = ['--X', r['X'], '--A', r['A']]
    return a + (['--ring', r['ring']] if r['ring'] else ['--Y', str(int(float(r['Y'])))])


def build(xyz, r, outdir, name, setname, extra=()):
    res = cli('build-xyz', xyz, '--out', outdir, '--name', name, '--nproc', RES['nproc'], '--mem', RES['mem'], *atomargs(r), *extra)
    for f in res['files']:
        prog = f.split('/')[0]
        if prog in ('gaussian', 'orca', 'qchem', 'qchem_eda', 'gaussian_cp', 'orca_cp') and not f.endswith('.pc'):
            jobs.append(dict(set=setname, name=name, program=prog, input=os.path.join(setname, f),
                             manifest=os.path.join(setname, 'manifests', name + '_manifest.json')))
    return res


def main(data, pkg):
    hal = list(csv.DictReader(open(os.path.join(data, 'hal59', 'hal59.csv'))))
    meta = {r['name']: r for r in hal}
    for r in hal:
        build(os.path.join(data, 'hal59', r['name'] + '.xyz'), r, os.path.join(pkg, 'A_HAL59'), r['name'], 'A_HAL59')
    for r in csv.DictReader(open(os.path.join(data, 'scans', 'scans.csv'))):
        sub = 'B_scan_distance' if r['scan'] == 'distance' else 'C_scan_angle'
        build(os.path.join(data, 'scans', r['file']), r, os.path.join(pkg, sub), r['file'][:-4], sub)
    for nm in SENS:
        r = meta[nm]; xyz = os.path.join(data, 'hal59', nm + '.xyz')
        for m in METHODS:
            build(xyz, r, os.path.join(pkg, 'D_method_basis'), f'{nm}__{m}__def2-TZVP', 'D_method_basis', ['--method', m, '--basis', 'def2-TZVP', '--noEda', '1'])
        for b in BASES:
            build(xyz, r, os.path.join(pkg, 'D_method_basis'), f'{nm}__wB97X-D__{b}', 'D_method_basis', ['--method', 'wB97X-D', '--basis', b, '--noEda', '1'])
        for b in ('def2-SVPD', 'def2-QZVP'):
            build(xyz, r, os.path.join(pkg, 'E_eda_settings'), f'{nm}__EDA__{b}', 'E_eda_settings', ['--programs', 'none', '--edaBasis', b])
        build(xyz, r, os.path.join(pkg, 'E_eda_settings'), f'{nm}__EDA__def2-TZVPD__noBSSE', 'E_eda_settings', ['--programs', 'none', '--noBSSE', '1'])
    with open(os.path.join(pkg, 'jobs_small.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, list(jobs[0].keys())); w.writeheader(); w.writerows(jobs)
    for k, v in sorted(Counter((j['set'], j['program']) for j in jobs).items()): print(k, v)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '../data', sys.argv[2] if len(sys.argv) > 2 else '../qm_package')
