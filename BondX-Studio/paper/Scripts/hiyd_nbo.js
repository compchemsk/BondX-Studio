#!/usr/bin/env node
/* hiyd_nbo.js - read the existing hIYD QM-cluster NBO outputs (Gaussian 16, NBO 3.1) with the BondX T3 parser.
   Usage: node hiyd_nbo.js <hIYD_for_BondX root> <out.json>
   For every US/C?/qm-halo/<RES>-nbo.log:
     - B.parseNBO reads the second-order table and the last geometry;
     - covalent bonds come from covalent radii (1.25 x sum);
     - fragment 1 is the molecule that holds the iodine, fragment 2 is the rest of the cluster;
     - every N/O/S atom of fragment 2 and every five- or six-membered ring of fragment 2 is a candidate;
     - for each candidate the T1 geometry (d, theta1) and the BondX gate value
       E2_core = sum of LP(Y) -> sigma*(C-I) (atoms) or pi(ring) -> sigma*(C-I) (rings) are stored;
     - the sum of LP -> sigma*(C-I) terms whose donor lies in fragment 1 (intramolecular) is stored separately.
   No new calculation is done. */
'use strict';
const fs = require('fs'), path = require('path');
const B = require('./engine.js');
const [root, outFile] = process.argv.slice(2);
const P = B.defaults(); const cut = P.geom.cutoffs.I;
const res = [];
for (const C of ['C1', 'C2', 'C3']) {
  const dir = path.join(root, 'US', C, 'qm-halo');
  for (const f of fs.readdirSync(dir).filter(x => x.endsWith('-nbo.log')).sort()) {
    const text = fs.readFileSync(path.join(dir, f), 'utf8');
    const nbo = B.parseNBO(text); const g = nbo.geometry; const n = g.length;
    const xyz = g.map(a => [a.x, a.y, a.z]);
    const bonded = Array.from({ length: n }, () => []);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const r = (B.COV[g[i].el] || 0.8) + (B.COV[g[j].el] || 0.8);
      if (B.dist(xyz[i], xyz[j]) < 1.25 * r) { bonded[i].push(j); bonded[j].push(i); }
    }
    const iI = g.findIndex(a => a.el === 'I'); if (iI < 0) throw new Error(f + ': no iodine');
    const nbC = bonded[iI].filter(j => g[j].el === 'C'); if (nbC.length !== 1) throw new Error(f + ': iodine has ' + nbC.length + ' carbon partners');
    const iC = nbC[0];
    const frag = new Array(n).fill(2); const st = [iI]; frag[iI] = 1;
    while (st.length) { const a = st.pop(); for (const b of bonded[a]) if (frag[b] !== 1) { frag[b] = 1; st.push(b); } }
    // candidates: atoms
    const xbs = []; const geo = [];
    g.forEach((a, k) => {
      if (frag[k] !== 2 || !['N', 'O', 'S'].includes(a.el)) return;
      const d = B.dist(xyz[iI], xyz[k]); const t1 = B.angle(xyz[iC], xyz[iI], xyz[k]);
      const nH = bonded[k].filter(j => g[j].el === 'H').length; const nHeavy = bonded[k].length - nH;
      xbs.push({ id: 'Y' + (k + 1), type: 'atom', X: iI + 1, A: iC + 1, Y: k + 1 });
      geo.push({ id: 'Y' + (k + 1), el: a.el, n: k + 1, d, t1, nH, nHeavy, t1ok: t1 >= 140 && t1 <= 180, dok: d <= cut[a.el] });
    });
    // candidates: 5/6-membered rings of fragment 2 (heavy atoms only)
    const heavy = k => frag[k] === 2 && g[k].el !== 'H'; const rings = new Set();
    const dfs = (start, cur, pathA) => {
      if (pathA.length > 6) return;
      for (const nb of bonded[cur]) {
        if (!heavy(nb)) continue;
        if (nb === start && pathA.length >= 5) { rings.add(pathA.slice().sort((a, b) => a - b).join(',')); continue; }
        if (pathA.includes(nb) || nb < start) continue;
        pathA.push(nb); dfs(start, nb, pathA); pathA.pop();
      }
    };
    for (let k = 0; k < n; k++) if (heavy(k)) dfs(k, k, [k]);
    for (const rk of rings) {
      const ids = rk.split(',').map(Number); const cen = [0, 1, 2].map(c => ids.reduce((s, i) => s + xyz[i][c], 0) / ids.length);
      const d = B.dist(xyz[iI], cen); const t1 = B.angle(xyz[iC], xyz[iI], cen);
      xbs.push({ id: 'R' + ids.map(i => i + 1).join('-'), type: 'pi', X: iI + 1, A: iC + 1, ring: ids.map(i => i + 1) });
      geo.push({ id: 'R' + ids.map(i => i + 1).join('-'), el: 'AR', n: null, size: ids.length, d, t1, t1ok: t1 >= 140 && t1 <= 180, dok: d <= cut.AR });
    }
    const gate = B.nboGate({ xbs, atoms: g.map((a, k) => ({ n: k + 1, el: a.el, frag: frag[k] })) }, nbo, { threshold: 0.5, includeRY: false });
    // all LP -> sigma*(C-I) entries, split by donor fragment
    let intra = 0, inter = 0; const lpRows = [];
    const hops = Array(n).fill(-1); { hops[iC] = 0; const qu = [iC]; while (qu.length) { const a = qu.shift(); for (const b of bonded[a]) if (hops[b] < 0) { hops[b] = hops[a] + 1; qu.push(b); } } }
    for (const e of nbo.entries) {
      const a = e.acceptor, d = e.donor;
      if (!(a.type === 'BD*' && a.atoms.length === 2 && a.atoms.some(t => t.n === iI + 1) && a.atoms.some(t => t.n === iC + 1))) continue;
      if (d.type !== 'LP') continue;
      const fr = frag[d.atoms[0].n - 1]; if (fr === 1) intra += e.e2; else inter += e.e2;
      lpRows.push({ donor: d.atoms[0].el + d.atoms[0].n, frag: fr, e2: e.e2, bondsToC: hops[d.atoms[0].n - 1], dToC: +B.dist(xyz[d.atoms[0].n - 1], xyz[iC]).toFixed(3) });
    }
    res.push({ C, residue: f.replace('-nbo.log', ''), terminated: nbo.terminated, nAtoms: n, nFrag1: frag.filter(x => x === 1).length, I: iI + 1, Cx: iC + 1,
      geo, gate: gate.map(x => ({ id: x.id, E2_core: x.E2_core, E2_total: x.E2_total, E2_RY: x.E2_RY, hbY: x.hbY, lpOnY: x.lpOnY, verdict: x.verdict, contrib: x.contrib })), lpIntra: intra, lpInter: inter, lpRows });
  }
}
fs.writeFileSync(outFile, JSON.stringify(res, null, 1));
process.stdout.write(JSON.stringify({ logs: res.length, terminated: res.filter(r => r.terminated).length }) + '\n');
