#!/usr/bin/env node
/* hiyd_bondx.js - run the BondX T1 detector on the hIYD umbrella-sampling windows.
   Usage: node hiyd_bondx.js <hIYD_for_BondX root> <out.jsonl>
   For every window US/C?/traj/C?_US-<w>.{pdb,dcd}:
     - the BondX reader removes the extra point (atom name EP) of the sigma-hole model and maps the DCD frames;
     - every DCD frame is tested with the BondX T1 criteria;
     - one JSON line per frame lists all T1 contacts of the iodine of residue MIA.
   Settings: BondX default cut-offs and theta1 window (140-180 deg), the other ligand (FMN, FRA) is an acceptor,
   water is absent from the stripped files, and the pi-normal limit is opened to 90 deg so that every ring
   contact is recorded. The stricter BondX filters (Lewis class, pi normal <= 60 deg) are applied later in Python
   from the stored fields. */
'use strict';
const fs = require('fs'), path = require('path');
const B = require('./engine.js');
const [root, outFile] = process.argv.slice(2);
const P = B.mergeParams(B.defaults(), { geom: { piNormalMax: 90, acceptorScope: { protein: true, nucleic: true, water: false, otherLigand: true }, excludeNonLewis: false, reportNearMiss: false } });
const out = fs.openSync(outFile, 'w');
const summary = [];
for (const C of ['C1', 'C2', 'C3']) {
  const dir = path.join(root, 'US', C, 'traj');
  const wins = fs.readdirSync(dir).filter(f => f.endsWith('.dcd')).map(f => f.replace(/\.dcd$/, ''));
  for (const stem of wins) {
    const w = parseFloat(stem.split('US-')[1]);
    const S = B.parsePDB(fs.readFileSync(path.join(dir, stem + '.pdb'), 'utf8'), { altloc: 'first' }); // the extra point EP is removed here (S.vsites)
    if (!S.vsites || !S.vsites.names.includes('MIA:EP')) throw new Error(stem + ': extra point of MIA not found');
    const buf = fs.readFileSync(path.join(dir, stem + '.dcd'));
    const D = B.parseDCD(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength), S.vsites.fileAtoms);
    const frames = B.dropVsites(S, D.frames);
    const prep = B.prepareDetection(S, P);
    prep.donors = prep.donors.filter(d => S.residues[d.res].resName === 'MIA'); // the unbinding ligand only (MIB sits in the other site)
    const donors = prep.donors.map(d => ({ X: S.atoms[d.X].name, A: S.atoms[d.A].name, res: S.residues[d.res].resName + S.residues[d.res].resSeq }));
    if (donors.length !== 1 || donors[0].X !== 'I1' || donors[0].A !== 'C6' || !donors[0].res.startsWith('MIA')) throw new Error(stem + ': unexpected donor ' + JSON.stringify(donors));
    frames.forEach((f, k) => {
      const hits = B.detectFrame(S, prep, f, D.boxes[k]).filter(h => h.status === 'XB').map(h => {
        B.decorate(S, h); const R = S.residues[h.accRes];
        return { res: R.resName, seq: R.resSeq, atom: h.accName, el: h.accEl, type: h.type, cls: h.accClass, lewis: h.lewis,
          d: +h.d.toFixed(4), t1: +h.theta1.toFixed(3), t2: h.theta2 == null ? null : +h.theta2.toFixed(3), G: +h.G.toFixed(4) };
      });
      fs.writeSync(out, JSON.stringify({ C, w, k, n: D.frames.length, hits }) + '\n');
    });
    summary.push({ C, w, frames: D.frames.length, atoms: S.atoms.length, donor: donors[0] });
  }
}
fs.closeSync(out);
process.stdout.write(JSON.stringify({ windows: summary.length, frames: summary.reduce((s, x) => s + x.frames, 0), bondxVersion: B.VERSION }) + '\n');
