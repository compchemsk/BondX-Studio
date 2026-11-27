#!/usr/bin/env node
/* 19_v3_field.js - MM field at the halogen without any QM job (review comments M11, M12, minor 12-14).
   Usage: node 19_v3_field.js ../data/protein ../results/v3
   For the six protein systems (same contacts and QM regions as the paper models) it computes the potential phi_X and the
   field component E_AX along A->X from the MM point charges for
     - embedding radii 10-40 A and the whole protein, with folded and with explicit MM hydrogens near X (mm.explicitH 6 A),
     - the RC and RCD boundary schemes,
     - ff14SB and CHARMM36 (number of MM atoms, of zero-charge atoms that are not written, and of written point charges). */
'use strict';
const fs = require('fs'), path = require('path');
const B = require('./engine.js');
const [PDB, OUT] = process.argv.slice(2);
const SYSTEMS = [['1us0_idd594', 'LDT:A:320'], ['2w0s', 'BVP:A:1207'], ['1p5e', 'TBS:A:301'], ['4agl', 'P84:A:400'], ['4gql', 'R47:A:306'], ['4dst', '9LI:A:201']];
const K = 14.399645;
const rows = [['system', 'contact', 'ff', 'boundary', 'radius', 'explicitH', 'nMM', 'nZero', 'nCharges', 'nExplicitH', 'net', 'phi_X', 'E_AX']];
for (const [sys, lig] of SYSTEMS) {
  const S = B.parsePDB(fs.readFileSync(path.join(PDB, sys + '.pdb'), 'utf8'));
  const base = { qm: { env: 'qmmm', keepAmide: 'off', orientPolarH: false, hbShell: 0 }, mm: { radius: 15, boundary: 'RC', explicitH: 0 }, ff: { name: 'ff14SB' } };
  const P0 = B.mergeParams(B.defaults(), base);
  const [rn, ch, sq] = lig.split(':');
  let xb = B.detectXB(S, P0).hits.filter(h => h.status === 'XB' && S.residues[h.donorRes].resName === rn && S.residues[h.donorRes].chain === ch && String(S.residues[h.donorRes].resSeq) === sq);
  xb = xb.map((h, k) => Object.assign({ cid: 'XB' + (k + 1) }, h));
  const variants = [];
  for (const r of [10, 15, 20, 25, 30, 40, 200]) for (const eh of [0, 6]) variants.push({ ff: 'ff14SB', boundary: 'RC', radius: r, explicitH: eh });
  variants.push({ ff: 'ff14SB', boundary: 'RCD', radius: 15, explicitH: 0 }, { ff: 'ff14SB', boundary: 'RCD', radius: 15, explicitH: 6 });
  variants.push({ ff: 'charmm36', boundary: 'RC', radius: 15, explicitH: 0 }, { ff: 'charmm36', boundary: 'RC', radius: 15, explicitH: 6 });
  for (const v of variants) {
    const P = B.mergeParams(P0, { mm: { radius: v.radius, boundary: v.boundary, explicitH: v.explicitH }, ff: { name: v.ff } });
    const ffa = B.assignForceField(S, P); const m = B.buildQMModel(S, 0, xb, { P, ffa });
    const pcs = B.embeddingCharges(S, 0, m.include, { P, ffa, qmResidueCharges: m.qmResCharge });
    for (const x of m.manifest.xbs) {
      const X = m.qm[x.X - 1], A = m.qm[x.A - 1]; const u = [X.x - A.x, X.y - A.y, X.z - A.z]; const un = Math.hypot(...u); u.forEach((c, i) => u[i] = c / un);
      let phi = 0, E = [0, 0, 0];
      for (const c of pcs.charges) { const r = [X.x - c.x, X.y - c.y, X.z - c.z]; const d = Math.hypot(...r); phi += K * c.q / d; for (let i = 0; i < 3; i++) E[i] += K * c.q * r[i] / (d * d * d); }
      rows.push([sys, x.label.replace(/ /g, ''), v.ff, v.boundary, v.radius === 200 ? 'full' : v.radius, v.explicitH, pcs.nMM, pcs.nZero, pcs.charges.length, pcs.nExplicitH, pcs.net.toFixed(3), phi.toFixed(4), (E[0] * u[0] + E[1] * u[1] + E[2] * u[2]).toFixed(4)]);
    }
  }
}
fs.writeFileSync(path.join(OUT, 'field_v3.csv'), rows.map(r => r.map(v => /[,"]/.test(String(v)) ? '"' + v + '"' : v).join(',')).join('\n') + '\n');
console.log(rows.length - 1, 'rows');
