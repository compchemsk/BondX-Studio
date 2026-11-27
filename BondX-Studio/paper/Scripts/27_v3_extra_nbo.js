#!/usr/bin/env node
/* 27_v3_extra_nbo.js - reads the Gaussian NBO outputs of the additional job sets with the BondX engine (version 3.0).
   G_hard_negatives : hIYD frames with poor-Lewis-base acceptors and Lewis-base controls (models built with version-3.0 defaults)
   J_4dst_amide     : 4DST Cl...N Gly75 with the amide of Gly75 intact (cluster gas, cluster PCM, QM/MM)
   K_hiyd_amide     : the 18 published hIYD backbone-N poses with the N-terminal H cap rebuilt as an acetyl group
   H_marginal_levels: written but not usable (every job stopped in its first link; duplicated route keywords)
   Usage: node 27_v3_extra_nbo.js <qm_extra folder> <out folder>
   Output: extra_nbo.csv (one row per contact and job) */
'use strict';
const fs = require('fs'), path = require('path');
const B = require('./engine.js');
const ROOT = process.argv[2] || '../qm_extra', OUT = process.argv[3] || '../results/v3';
// CSV with quoted fields (the role column contains commas)
const csvLine = l => { const out = []; let cur = '', inq = false; for (let i = 0; i < l.length; i++) { const c = l[i];
  if (inq) { if (c === '"' && l[i + 1] === '"') { cur += '"'; i++; } else if (c === '"') inq = false; else cur += c; }
  else if (c === '"') inq = true; else if (c === ',') { out.push(cur); cur = ''; } else cur += c; } out.push(cur); return out; };
const jobs = {}; { const L = fs.readFileSync(path.join(ROOT, 'jobs_extra.csv'), 'utf8').split(/\r?\n/).filter(Boolean); const h = csvLine(L[0]);
  L.slice(1).forEach(l => { const v = csvLine(l); const o = {}; h.forEach((c, i) => o[c] = v[i]); jobs[o.name] = o; }); }
const rows = [['set', 'name', 'role', 'contact', 'frame', 'd', 'theta1', 'G', 'lewis', 'terminated', 'nbo_complete', 'geometry_ok', 'E2_valence', 'E2_RY', 'verdict', 'hbY', 'hbFrag2', 'flags', 'lpOnY', 'intra_O_sigmaCI', 'E_pub_gate']];
const f3 = v => (v === null || v === undefined || Number.isNaN(v)) ? '' : (+v).toFixed(3);
const q = s => '"' + String(s).replace(/"/g, '""') + '"';
for (const set of ['G_hard_negatives', 'H_marginal_levels', 'J_4dst_amide', 'K_hiyd_amide']) {
  const dir = path.join(ROOT, set, 'gaussian'); if (!fs.existsSync(dir)) continue;
  for (const f of fs.readdirSync(dir).filter(x => x.endsWith('.gjf')).sort()) {
    const stem = f.replace(/\.gjf$/, ''); const J = jobs[stem] || {}; const log = path.join(dir, stem + '.log');
    const base = [set, stem, q(J.role || ''), q(J.contact || ''), J.frame || '', J.d || '', J.theta1 || '', J.G || '', J.lewis || ''];
    if (!fs.existsSync(log)) { rows.push(base.concat(['missing'])); continue; }
    const txt = fs.readFileSync(log, 'utf8'); const nbo = B.parseNBO(txt);
    const complete = /Natural Bond Orbitals \(Summary\)/.test(txt);
    if (!nbo.entries.length) { rows.push(base.concat([nbo.terminated, complete])); continue; }
    const mf = path.join(ROOT, set, 'manifests', stem + '_manifest.json');
    let man, geoOK = '';
    if (fs.existsSync(mf)) { man = JSON.parse(fs.readFileSync(mf, 'utf8')); geoOK = B.verifyGeometry(man, nbo.geometry).ok; }
    else { // set K: the acceptor is atom 1 (backbone N); the donor is the C-I bond of I-Tyr; the I-Tyr atoms form fragment 1
      const geo = nbo.geometry; const iI = geo.findIndex(a => a.el === 'I') + 1; const pI = geo[iI - 1];
      let iC = -1, dmin = 9; geo.forEach((a, k) => { if (a.el !== 'C') return; const d = Math.hypot(a.x - pI.x, a.y - pI.y, a.z - pI.z); if (d < dmin) { dmin = d; iC = k + 1; } });
      // fragment 1 = atoms connected to the iodine through bonds (covalent radii), so that its own lone pairs are not fragment 2
      const R = { H: 0.31, C: 0.76, N: 0.71, O: 0.66, S: 1.05, I: 1.39 }; const nA = geo.length; const frag = new Array(nA).fill(2); const st = [iI - 1]; frag[iI - 1] = 1;
      while (st.length) { const a = st.pop(); for (let b = 0; b < nA; b++) { if (frag[b] === 1) continue; const d = Math.hypot(geo[a].x - geo[b].x, geo[a].y - geo[b].y, geo[a].z - geo[b].z); if (d < (R[geo[a].el] || 0.8) + (R[geo[b].el] || 0.8) + 0.45) { frag[b] = 1; st.push(b); } } }
      man = { atoms: geo.map((a, k) => ({ n: k + 1, el: a.el, frag: frag[k] })), fragCharge: {}, xbs: [{ id: 'XB1', type: 'atom', X: iI, A: iC, Y: 1, label: 'I...N backbone (acetyl cap)', lewis: 'weak' }] };
      geoOK = 'no manifest';
    }
    const intra = nbo.entries.filter(e => e.donor.type === 'LP' && e.donor.atoms[0].el === 'O' && man.atoms[e.donor.atoms[0].n - 1].frag === 1 && e.acceptor.type === 'BD*' && e.acceptor.atoms.some(t => t.el === 'I')).reduce((m, e) => Math.max(m, e.e2), 0);
    B.nboGate(man, nbo).forEach(g => rows.push(base.slice(0, 3).concat([q(J.contact || g.label)], base.slice(4), [nbo.terminated, complete, geoOK, f3(g.E2_valence), f3(g.E2_RY), g.verdict, f3(g.hbY), f3(g.hbFrag2), g.flags.join('|'), g.lpOnY, f3(intra), J.E_pub_gate || ''])));
  }
}
fs.mkdirSync(OUT, { recursive: true }); fs.writeFileSync(path.join(OUT, 'extra_nbo.csv'), rows.map(r => r.join(',')).join('\n') + '\n');
console.log(`extra_nbo.csv: ${rows.length - 1} rows`);
