#!/usr/bin/env python3
"""14_ligand_check.py
Checks the ligand protonation of BondX (one heavy-atom structure for hydrogens and charge) against the PDB
FORMUL record of every ligand residue in the protein test set (56 files).
Expected H = FORMUL H + (model charge - FORMUL charge) - (covalent links to other residues).
Writes results/ligand_check.csv and prints a summary.
Usage: python3 14_ligand_check.py ..
"""
import os, sys, csv, json
from bxutil import cli

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
PDB = os.path.join(ROOT, 'data', 'protein'); RES = os.path.join(ROOT, 'results')
rows = []
for fn in sorted(f for f in os.listdir(PDB) if f.endswith('.pdb')):
    d = cli('ligand-check', os.path.join(PDB, fn))
    for r in d['ligands']: rows.append(dict(pdb=fn[:-4], **r))
with open(os.path.join(RES, 'ligand_check.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)
metal = [r for r in rows if r['metal']]; nof = [r for r in rows if not r['metal'] and not r['formul']]
inc = [r for r in rows if r['formul'] and not r['complete']]; comp = [r for r in rows if r['formul'] and r['complete']]
ok = [r for r in comp if r['hOK']]
summary = dict(files=len({r['pdb'] for r in rows}), ligand_residues=len(rows), metal=len(metal), no_formul=len(nof), incomplete=len(inc),
               complete=len(comp), H_matches=len(ok), mismatches=sorted({f"{r['pdb']} {r['resName']}" for r in comp if not r['hOK']}))
json.dump(summary, open(os.path.join(RES, 'ligand_check_summary.json'), 'w'), indent=1)
print(json.dumps(summary, indent=1))
