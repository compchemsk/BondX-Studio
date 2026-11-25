#!/usr/bin/env python3
"""Independent re-implementation (numpy only, no BondX code) of the published halogen-bond criteria
(halo-bond.py of the hIYD study) on the stripped DCD files. Compares every frame with results/hiyd/bondx_frames.jsonl."""
import os, re, json, struct, glob, sys
import numpy as np
import sys
ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
D = sys.argv[2] if len(sys.argv) > 2 else '/path/to/hIYD_for_BondX/US'
P = ROOT + '/results/hiyd/bondx_frames.jsonl'
CUT = {'N': 3.97, 'O': 3.82, 'S': 4.15, 'aromatic': 4.00}
RING = {'TYR': [['CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ']], 'PHE': [['CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ']],
        'TRP': [['CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'CE2'], ['CE2', 'CD2', 'CG', 'CD1', 'NE1']]}

def read_dcd(fn):
    b = open(fn, 'rb').read(); o = 0
    def rec():
        nonlocal o
        n = struct.unpack_from('<i', b, o)[0]; d = b[o + 4:o + 4 + n]; o += 8 + n; return d
    h = rec(); nset = struct.unpack_from('<i', h, 4)[0]; charmm = struct.unpack_from('<i', h, 80)[0] != 0; uc = charmm and struct.unpack_from('<i', h, 44)[0] == 1
    rec(); natom = struct.unpack('<i', rec())[0]; fr = []
    while o < len(b):
        if uc: rec()
        x = np.frombuffer(rec(), '<f4'); y = np.frombuffer(rec(), '<f4'); z = np.frombuffer(rec(), '<f4')
        fr.append(np.stack([x, y, z], 1).astype(np.float64))
    return fr

tot = mism = nf = 0; per = {}
bx = {}
for l in open(P):
    r = json.loads(l); bx[(r['C'], r['w'], r['k'])] = r['hits']
for C in ['C1', 'C2', 'C3']:
    for pdb in sorted(glob.glob(f'{D}/{C}/traj/*.pdb')):
        w = float(pdb.split('US-')[1][:-4]); atoms = []
        for l in open(pdb):
            if l.startswith(('ATOM', 'HETATM')): atoms.append((l[12:16].strip(), l[17:21].strip(), int(l[22:26])))
        iI = [i for i, a in enumerate(atoms) if a[1] == 'MIA' and a[0] == 'I1'][0]; iC = [i for i, a in enumerate(atoms) if a[1] == 'MIA' and a[0] == 'C6'][0]
        acc = []
        res = {}
        for i, a in enumerate(atoms): res.setdefault((a[1], a[2]), {})[a[0]] = i
        for (rn, rs), names in res.items():
            if rn in ('MIA',): continue
            for grp in RING.get(rn, []):
                ids = [names[n] for n in grp if n in names]
                if ids: acc.append((ids, 'aromatic', rn, rs, 'ring'))
            for n, i in names.items():
                if n == 'EP': continue
                if n[0] in 'NOS': acc.append((i, n[0], rn, rs, n))
        for k, X in enumerate(read_dcd(pdb[:-4] + '.dcd')):
            nf += 1; hits = set()
            pI, pC = X[iI], X[iC]
            for a in acc:
                p = X[a[0]].mean(0) if a[1] == 'aromatic' else X[a[0]]
                d = np.linalg.norm(p - pI)
                if d <= CUT[a[1]]:
                    v1, v2 = pC - pI, p - pI; ang = np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1)))
                    if 140 <= ang <= 180: hits.add((a[2], a[3], a[1], a[4]))
            b = set()
            for h in bx[(C, w, k)]:
                rn = h['res']; b.add((rn, h['seq'], 'aromatic' if h['type'] == 'pi' else h['el'], 'ring' if h['type'] == 'pi' else h['atom']))
            # the stripped PDB may name HIS as HIE/HID: names are compared as written in the file on both sides
            tot += len(hits)
            if hits != b: mism += 1; per[(C, w, k)] = (sorted(hits ^ b))
print(json.dumps({'frames': nf, 'contacts_independent': tot, 'frames_differing': mism, 'examples': [str(x) for x in list(per.items())[:5]]}, indent=1))
