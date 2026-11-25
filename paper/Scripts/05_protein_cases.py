#!/usr/bin/env python3
"""05_protein_cases.py
Protein-ligand test set. BondX T1 detection and T3a QM models (cluster and QM/MM) through bondx_cli.js.
Structures: PDB files distributed with the PLIP test suite (literature-validated halogen bonds) and 1US0 + IDD594.
For every model variant the script writes Gaussian 16, ORCA and Q-Chem NBO inputs (and a Q-Chem ALMO-EDA input for
cluster models), and it records model size, charges, electron count and MM-embedding statistics.
It also computes, from the MM point charges alone, the electrostatic potential (ESP) at the halogen X and the field
component along the A->X bond axis. These two numbers show how much the environment changes between force fields,
boundary schemes and embedding radii before any QM job is run.
Output: qm_package/F_protein/..., results/protein_models.csv, results/protein_t1.csv
"""
import os, sys, csv, json
import numpy as np
from bxutil import cli

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
PDB = os.path.join(ROOT, 'data', 'protein')
PKG = os.path.join(ROOT, 'qm_package', 'F_protein')
SYSTEMS = [('1us0_idd594', 'LDT:A:320'), ('2w0s', 'BVP:A:1207'), ('1p5e', 'TBS:A:301'), ('4agl', 'P84:A:400'), ('4gql', 'R47:A:306'), ('4dst', '9LI:A:201')]
REPLICAS = [('2w0s', 'BVP:B:1207'), ('1p5e', 'TBS:C:302'), ('4agl', 'P84:B:400')]
CL = {'qm': {'env': 'cluster', 'pcm': True, 'eps': 4}}
QMMM = lambda **k: {'qm': {'env': 'qmmm'}, 'mm': dict({'radius': 15, 'boundary': 'RC'}, **k.get('mm', {})), 'ff': dict({'name': 'ff14SB'}, **k.get('ff', {}))}
VARIANTS = [
    ('cluster_pcm4', CL, True),
    ('cluster_gas', {'qm': {'env': 'cluster', 'pcm': False}}, False),
    ('qmmm_ff14SB', QMMM(), False),
    ('qmmm_charmm36', QMMM(ff={'name': 'charmm36'}), False),
    ('qmmm_Z1', QMMM(mm={'boundary': 'Z1'}), False),
    ('qmmm_Z2', QMMM(mm={'boundary': 'Z2'}), False),
    ('qmmm_Z3', QMMM(mm={'boundary': 'Z3'}), False),
    ('qmmm_r10', QMMM(mm={'radius': 10}), False),
    ('qmmm_r20', QMMM(mm={'radius': 20}), False),
    ('qmmm_r25', QMMM(mm={'radius': 25}), False),
    ('qmmm_full', QMMM(mm={'radius': 200}), False),
    ('shell3', {'qm': {'env': 'cluster', 'pcm': True, 'eps': 4, 'shell': 3}}, False),
    ('shell4', {'qm': {'env': 'cluster', 'pcm': True, 'eps': 4, 'shell': 4}}, False),
    ('shell5', {'qm': {'env': 'cluster', 'pcm': True, 'eps': 4, 'shell': 5}}, False),
]
EXTRA = [('1us0_idd594', 'LDT:A:320', 'qmmm_HIS110_HID', QMMM(ff={'overrides': {'A:110': 'HID'}})),
         ('1us0_idd594', 'LDT:A:320', 'qmmm_HIS110_HIP', QMMM(ff={'overrides': {'A:110': 'HIP'}}))]
# ff19SB uses the ff14SB protein charges (it differs in backbone torsions/CMAP), so its point charges are identical.
# The script checks this for every system (results/ff19SB_check.csv) instead of writing duplicate QM jobs.
MAX_ATOMS = 260   # larger shell models are not written (too costly for NBO at def2-TZVP)
jobs, rows = [], []


def esp_field(manifest, csvfile):
    """ESP (V) at X and field component along A->X (V/A) from the MM point charges."""
    if not os.path.exists(csvfile): return None, None
    q = np.loadtxt(csvfile, delimiter=',', skiprows=1, usecols=(0, 1, 2, 3), ndmin=2)
    out = []
    for xb in manifest['xbs']:
        X = np.array([manifest['atoms'][xb['X'] - 1][k] for k in 'xyz']); A = np.array([manifest['atoms'][xb['A'] - 1][k] for k in 'xyz'])
        r = X - q[:, :3]; d = np.linalg.norm(r, axis=1); k = 14.399645  # e/(4 pi eps0) in V*A
        phi = k * np.sum(q[:, 3] / d); E = k * np.sum((q[:, 3] / d ** 3)[:, None] * r, axis=0)
        u = (X - A) / np.linalg.norm(X - A)
        out.append((phi, float(E @ u)))
    return out[0]


def record(name, sysname, lig, variant, res):
    e = res.get('embedding') or {}
    mf = json.load(open(os.path.join(PKG, 'manifests', name + '_manifest.json')))
    ef = esp_field(mf, os.path.join(PKG, 'mm_charges', name + '_mm_charges.csv')) if e else (None, None)
    rows.append(dict(model=name, system=sysname, ligand=lig, variant=variant, n_atoms=res['nAtoms'], n_heavy=res['nHeavy'], added_H=res['addedH'], link_H=res['linkH'],
                     charge_donor=res['fragCharge']['1'], charge_acceptor=res['fragCharge']['2'], total_charge=res['totalCharge'], multiplicity=res['multiplicity'], electrons=res['electrons'],
                     n_xb=len(res['xbs']), xbs=' ; '.join(x['label'] for x in res['xbs']), mm_n=e.get('n', ''), mm_net=round(e['net'], 4) if e else '', mm_M1=e.get('M1', ''),
                     mm_cut=e.get('cut', ''), boundary=e.get('boundary', ''), radius=e.get('radius', ''), esp_X_V=round(ef[0], 4) if ef[0] is not None else '', field_AX_VperA=round(ef[1], 4) if ef[1] is not None else '',
                     warnings=' | '.join(res['warnings'] + (e.get('notes') or []))))
    for f in res['files']:
        prog = f.split('/')[0]
        if prog in ('gaussian', 'orca', 'qchem', 'qchem_eda', 'gaussian_cp', 'orca_cp') and not f.endswith('.pc'):
            jobs.append(dict(set='F_protein', name=name, program=prog, input=os.path.join('F_protein', f), manifest=os.path.join('F_protein', 'manifests', name + '_manifest.json')))


def build(sysname, lig, variant, P, eda):
    name = f'{sysname}__{lig.split(":")[0]}{lig.split(":")[1]}__{variant}'
    args = ['build-pdb', os.path.join(PDB, sysname + '.pdb'), '--ligand', lig, '--out', PKG, '--name', name, '--P', json.dumps(P), '--nproc', 40, '--mem', 50000, '--cpRobust', '1']
    if not eda: args += ['--noEda', '1']
    probe = cli(*(args[:4] + ['--out', '/tmp/bx_probe', '--name', 'probe', '--P', json.dumps(P), '--programs', 'none', '--noEda', '1']))
    if probe.get('error'): print('skip', name, probe['error']); return
    if probe['nAtoms'] > MAX_ATOMS: print('skip', name, probe['nAtoms'], 'atoms'); rows.append(dict(model=name, system=sysname, ligand=lig, variant=variant, n_atoms=probe['nAtoms'], warnings='not written: more than %d atoms' % MAX_ATOMS)); return
    record(name, sysname, lig, variant, cli(*args))


def main():
    # T1 on every structure of the test set
    t1 = []
    for fn in sorted(os.listdir(PDB)):
        if not fn.endswith('.pdb'): continue
        for incF in (False, True):
            d = cli('detect-pdb', os.path.join(PDB, fn), '--P', json.dumps({'geom': {'reportNearMiss': True, 'includeF': incF}}))
            for h in d['hits']:
                t1.append(dict(pdb=fn[:-4], includeF=incF, atoms=d['atoms'], ms=round(d['msParse'] + d['msDetect'], 1), **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in h.items()}))
            if not d['hits']:
                t1.append(dict(pdb=fn[:-4], includeF=incF, atoms=d['atoms'], ms=round(d['msParse'] + d['msDetect'], 1), donor='; '.join(l['resName'] + ':' + l['chain'] + ':' + str(l['resSeq']) for l in d['halogenatedLigands']) or 'no halogenated ligand', status='none'))
    keys = sorted({k for r in t1 for k in r}, key=lambda k: (list(t1[0].keys()).index(k) if k in t1[0] else 99, k))
    with open(os.path.join(ROOT, 'results', 'protein_t1.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(t1)
    # QM models
    for s, lig in SYSTEMS:
        for v, P, eda in VARIANTS: build(s, lig, v, P, eda)
    for s, lig in REPLICAS: build(s, lig, 'cluster_pcm4', CL, True)
    chk = []
    for s, lig in SYSTEMS:
        a = cli('build-pdb', os.path.join(PDB, s + '.pdb'), '--ligand', lig, '--out', '/tmp/bx_ff19', '--name', 'a', '--P', json.dumps(QMMM()), '--programs', 'none', '--noEda', '1')
        b = cli('build-pdb', os.path.join(PDB, s + '.pdb'), '--ligand', lig, '--out', '/tmp/bx_ff19', '--name', 'b', '--P', json.dumps(QMMM(ff={'name': 'ff19SB'})), '--programs', 'none', '--noEda', '1')
        qa = open('/tmp/bx_ff19/mm_charges/a_mm_charges.csv').read().split('\n')[1:]; qb = open('/tmp/bx_ff19/mm_charges/b_mm_charges.csv').read().split('\n')[1:]
        same = [x.split(',')[:4] for x in qa] == [x.split(',')[:4] for x in qb]   # coordinates and charges; atom types differ by design
        chk.append(dict(system=s, n_ff14SB=a['embedding']['n'], n_ff19SB=b['embedding']['n'], identical_charges=same))
    with open(os.path.join(ROOT, 'results', 'ff19SB_check.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, list(chk[0].keys())); w.writeheader(); w.writerows(chk)
    print(chk)
    for s, lig, v, P in EXTRA: build(s, lig, v, P, False)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(ROOT, 'results', 'protein_models.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    with open(os.path.join(ROOT, 'qm_package', 'jobs_protein.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, list(jobs[0].keys())); w.writeheader(); w.writerows(jobs)
    print(len(rows), 'models', len(jobs), 'jobs')


if __name__ == '__main__':
    main()
