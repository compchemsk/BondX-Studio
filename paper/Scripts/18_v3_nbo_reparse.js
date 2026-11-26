#!/usr/bin/env node
/* 18_v3_nbo_reparse.js - re-reads every Gaussian NBO output with the BondX 3.0 parser (no new QM job).
   Usage: node 18_v3_nbo_reparse.js ../qm_package ../results/v3
   Writes
     nbo_v3.csv        one row per job and contact: valence SigmaE(2) (decision value), Rydberg part, verdict (grey zone),
                       LP(X)->sigma*(N-H/O-H) donation to the acceptor and to fragment 2, LP on the acceptor, NPA charge transfer
     pairs_v3.csv      every halogen...N/O/S pair (d <= 4.5 A) between fragment 1 and fragment 2 in the default protein models
                       (cluster_pcm4), with its T1 status and its valence SigmaE(2): the in-cluster negative set
     amideN_v3.csv     every N atom of an intact amide (N bonded to a C that carries a terminal O) in the protein models,
                       and whether NBO gives it a lone pair (review comment M9) */
'use strict';
const fs = require('fs'), path = require('path');
const B = require('./engine.js');
const [PKG, OUT] = process.argv.slice(2);
const P0 = B.defaults();
const jobs = [];
for (const f of ['jobs_small.csv', 'jobs_protein.csv']) {
  const L = fs.readFileSync(path.join(PKG, f), 'utf8').trim().split(/\r?\n/); const h = L[0].split(',');
  L.slice(1).forEach(l => { const t = l.split(','); const o = {}; h.forEach((k, i) => o[k] = t[i]); if (o.program === 'gaussian') jobs.push(o); });
}
const csv = rows => rows.map(r => r.map(v => v === null || v === undefined ? '' : (typeof v === 'string' && /[,"]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : (typeof v === 'number' ? +v.toFixed(4) : v))).join(',')).join('\n') + '\n';
const nboRows = [['set', 'name', 'xb', 'label', 'type', 'E2_valence', 'E2_RY', 'E2_with_RY', 'verdict', 'verdict_with_RY', 'hbY', 'hbFrag2', 'lpOnY', 'ct_NPA', 'geometry_ok']];
const pairRows = [['name', 'X', 'Y', 'Yname', 'Yres', 'Yel', 'd', 'theta1', 'cutoff', 'pass_dist', 'pass_angle', 'pass_T1', 'is_target', 'lewis_poor', 'E2_valence', 'E2_RY', 'hbY']];
const amideRows = [['name', 'N', 'src', 'lp_N', 'max_LP_N_to_sigmaCX']];
const el = s => s.toUpperCase();
for (const j of jobs) {
  const out = path.join(PKG, j.input.replace(/\.gjf$/, '.log')); if (!fs.existsSync(out)) continue;
  const txt = fs.readFileSync(out, 'utf8'); const nbo = B.parseNBO(txt); if (!nbo.entries.length) continue;
  const man = JSON.parse(fs.readFileSync(path.join(PKG, j.manifest), 'utf8'));
  const ver = B.verifyGeometry(man, nbo.geometry); const fq = B.fragmentNPA(man, nbo);
  const g0 = B.nboGate(man, nbo); const g1 = B.nboGate(man, nbo, { includeRY: true });
  g0.forEach((g, k) => nboRows.push([j.set, j.name, g.id, g.label, man.xbs[k].type, g.E2_valence, g.E2_RY, g1[k].E2_core, g.verdict, g1[k].verdict, g.hbY, g.hbFrag2, g.lpOnY, fq ? fq.ct : null, ver.ok]));
  if (j.set !== 'F_protein') continue;
  const A = man.atoms; const pos = n => [A[n - 1].x, A[n - 1].y, A[n - 1].z]; const d3 = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
  const ang = (a, b, c) => { const u = [a[0] - b[0], a[1] - b[1], a[2] - b[2]], v = [c[0] - b[0], c[1] - b[1], c[2] - b[2]]; return Math.acos(Math.max(-1, Math.min(1, (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (Math.hypot(...u) * Math.hypot(...v))))) * 180 / Math.PI; };
  // connectivity from distances (covalent radii + 0.45 A)
  const nbr = A.map(() => []); for (let i = 0; i < A.length; i++) for (let k = i + 1; k < A.length; k++) { const r = (B.COV[A[i].el] || 0.8) + (B.COV[A[k].el] || 0.8) + 0.45; if (d3(pos(i + 1), pos(k + 1)) < r) { nbr[i].push(k + 1); nbr[k].push(i + 1); } }
  // amide N atoms (any fragment)
  for (let i = 0; i < A.length; i++) {
    if (A[i].el !== 'N') continue; const n = i + 1;
    const amideC = nbr[i].find(c => A[c - 1].el === 'C' && nbr[c - 1].some(o => A[o - 1].el === 'O' && nbr[o - 1].length === 1));
    if (!amideC) continue;
    const lp = nbo.entries.some(e => e.donor.type === 'LP' && e.donor.atoms.length === 1 && e.donor.atoms[0].n === n);
    const toX = nbo.entries.filter(e => e.donor.type === 'LP' && e.donor.atoms[0].n === n && e.acceptor.type === 'BD*' && e.acceptor.atoms.some(t => ['CL', 'BR', 'I'].includes(el(t.el)))).reduce((m, e) => Math.max(m, e.e2), 0);
    amideRows.push([j.name, n, A[i].src, lp, toX]);
  }
  if (!/__cluster_pcm4$/.test(j.name)) continue;
  // every halogen (fragment 1, bonded to C) against every N/O/S of fragment 2
  for (let i = 0; i < A.length; i++) {
    if (!['CL', 'BR', 'I'].includes(A[i].el) || A[i].frag !== 1) continue; const X = i + 1;
    const Cx = nbr[i].find(c => A[c - 1].el === 'C'); if (!Cx) continue;
    for (let k = 0; k < A.length; k++) {
      if (A[k].frag === 1 || !['N', 'O', 'S'].includes(A[k].el)) continue; const Y = k + 1; const d = d3(pos(X), pos(Y)); if (d > 4.5) continue;
      if (A[k].src.startsWith('H(') || A[k].src.startsWith('link')) continue;
      const th = ang(pos(Cx), pos(X), pos(Y)); const cut = P0.geom.cutoffs[A[i].el][A[k].el];
      const passD = d <= cut, passA = th >= 140 && th <= 180;
      const pm = { xbs: [{ id: 'p', type: 'atom', X, A: Cx, Y, label: '' }], atoms: A, fragCharge: man.fragCharge };
      const g = B.nboGate(pm, nbo)[0];
      const isT = man.xbs.some(x => x.X === X && x.Y === Y);
      const nm = A[k].src.split('@')[1] || ''; const rn = (A[k].src.match(/^([A-Z0-9]+?)\d+/) || [])[1] || '';
      const poor = (nm === 'N') || (rn === 'LYS' && nm === 'NZ') || (rn === 'ARG' && ['NE', 'NH1', 'NH2'].includes(nm)) || (rn === 'TRP' && nm === 'NE1') || (rn === 'ASN' && nm === 'ND2') || (rn === 'GLN' && nm === 'NE2');
      pairRows.push([j.name, A[i].src, A[k].src, nm, A[k].src.split('@')[0], A[k].el, d, th, cut, passD, passA, passD && passA, isT, poor, g.E2_valence, g.E2_RY, g.hbY]);
    }
  }
}
fs.writeFileSync(path.join(OUT, 'nbo_v3.csv'), csv(nboRows)); fs.writeFileSync(path.join(OUT, 'pairs_v3.csv'), csv(pairRows)); fs.writeFileSync(path.join(OUT, 'amideN_v3.csv'), csv(amideRows));
console.log(JSON.stringify({ nbo: nboRows.length - 1, pairs: pairRows.length - 1, amideN: amideRows.length - 1, bondx: B.VERSION }));
