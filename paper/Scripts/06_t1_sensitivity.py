#!/usr/bin/env python3
"""06_t1_sensitivity.py
Sensitivity of BondX tier T1 to its two main criteria.
  theta1Min  : 130, 135, ..., 160 deg (default 140)
  cut-offs   : all X...Y distance cut-offs shifted by -0.2 ... +0.2 A (default 0)
For each setting: number of halogen bonds in the protein test set, recovery of the four literature-validated
contacts (2W0S Asn65, 1P5E Ile10 and Leu83, 4AGL Leu145; PLIP test suite), and the number of HAL59 halogen-bond
dimers that pass.
Output: results/t1_sensitivity.csv
"""
import os, sys, csv, json, copy
from bxutil import cli

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
PDB = os.path.join(ROOT, 'data', 'protein')
LIT = [('2w0s', 'BVP', 'A', 65), ('1p5e', 'TBS', 'A', 10), ('1p5e', 'TBS', 'A', 83), ('4agl', 'P84', 'A', 145)]
BASE_CUT = {'I': {'N': 3.97, 'O': 3.82, 'S': 4.15, 'AR': 4.00}, 'BR': {'N': 3.844, 'O': 3.6812, 'S': 4.02, 'AR': 3.87},
            'CL': {'N': 3.772, 'O': 3.6092, 'S': 3.948, 'AR': 3.798}, 'F': {'N': 3.32, 'O': 3.29, 'S': 3.57, 'AR': 3.40}}


def main():
    pdbs = sorted(f for f in os.listdir(PDB) if f.endswith('.pdb'))
    hal = [r for r in csv.DictReader(open(os.path.join(ROOT, 'results', 't1_hal59.csv'))) if r['role'] == 'halogen bond' and r['covered'] == 'True']
    rows = []
    for th in (130, 135, 140, 145, 150, 155, 160):
        for shift in (-0.2, -0.1, 0.0, 0.1, 0.2):
            cut = {x: {y: round(v + shift, 4) for y, v in d.items()} for x, d in BASE_CUT.items()}
            P = {'geom': {'theta1Min': th, 'cutoffs': cut}}
            n_xb, found = 0, set()
            for f in pdbs:
                d = cli('detect-pdb', os.path.join(PDB, f), '--P', json.dumps(P))
                xb = [h for h in d['hits'] if h['status'] == 'XB']; n_xb += len(xb)
                for pdb, lig, ch, res in LIT:
                    if f[:-4] == pdb and any(h['donorResName'] == lig and h['accChain'] == ch and h['accResSeq'] == res for h in xb): found.add((pdb, res))
            n_hal = 0
            for r in hal:
                x, y = r['halogen'], ('AR' if r['acc_el'] == 'pi' else r['acc_el'])
                if float(r['d']) <= cut[x][y] and float(r['theta1']) >= th: n_hal += 1
            rows.append(dict(theta1Min=th, cutoff_shift=shift, protein_XB=n_xb, literature_found=len(found), literature_total=len(LIT), hal59_pass=n_hal, hal59_covered=len(hal)))
            print(rows[-1])
    with open(os.path.join(ROOT, 'results', 't1_sensitivity.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, list(rows[0].keys())); w.writeheader(); w.writerows(rows)


if __name__ == '__main__':
    main()
