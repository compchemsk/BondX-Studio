#!/usr/bin/env node
/* 20_v3_hiyd_t2.js - tier T2 on the 108 hIYD umbrella-sampling windows (review comment M7).
   Usage: node 20_v3_hiyd_t2.js <hIYD_for_BondX root> ../results/v3
   Every window: BondX analyzeTrajectory with the default T1 criteria, FMN as an acceptor, no occupancy cut-off
   (dyn.occCut = 0, so that rare pocket contacts stay visible), and the real time between kept frames
   (stride from collect.log x 20 ps saving interval). Only the unbinding I-Tyr (MIA) is a donor.
   The windows are biased (umbrella restraint on xi); occupancies are per-window counts without reweighting. */
'use strict';
const fs = require('fs'), path = require('path');
const B = require('./engine.js');
const [root, OUT] = process.argv.slice(2);
const stride = {}; fs.readFileSync(path.join(root, 'collect.log'), 'utf8').split('\n').forEach(l => { const m = /\s+(C\d)_US-([\d.]+): (\d+) frames, stride (\d+)/.exec(l); if (m) stride[m[1] + '|' + (+m[2]).toFixed(1)] = +m[4]; });
const rows = [['pathway', 'xi', 'stride', 'dt_ps', 'frames', 'acceptor', 'atom', 'lewis', 'lewisWhy', 'occupancy', 'nEvents', 'tau_ps', 'tau_max_ps', 'D', 'bestFrame', 'bestD', 'bestT1', 'G_best']];
const PATH = { C1: 'A', C2: 'B', C3: 'C' };
for (const C of ['C1', 'C2', 'C3']) {
  const dir = path.join(root, 'US', C, 'traj');
  for (const f of fs.readdirSync(dir).filter(x => x.endsWith('.dcd')).sort()) {
    const stem = f.replace(/\.dcd$/, ''); const w = (+stem.split('US-')[1]).toFixed(1); const s = stride[C + '|' + w]; if (!s) throw new Error('no stride for ' + stem);
    const S = B.parsePDB(fs.readFileSync(path.join(dir, stem + '.pdb'), 'utf8'));
    const buf = fs.readFileSync(path.join(dir, f)); const D = B.parseDCD(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength), S.vsites.fileAtoms);
    const frames = B.dropVsites(S, D.frames);
    const P = B.mergeParams(B.defaults(), { geom: { acceptorScope: { protein: true, nucleic: true, water: false, otherLigand: true } }, dyn: { dt: s * 20, occCut: 0 } });
    const r = B.analyzeTrajectory(S, frames, D.boxes, P);
    r.results.filter(x => S.residues[x.donorRes].resName === 'MIA').forEach(x => {
      const gs = B.geometricScore(x.halogen, x.accEl === 'π' ? 'AR' : x.accEl, x.bestD, x.bestT1, x.type === 'pi' ? null : x.bestT2, P, x.type === 'pi', x.type === 'pi' ? x.bestT2 : null);
      rows.push([PATH[C], w, s, s * 20, frames.length, x.accLabel, x.accName, x.lewis, x.lewisWhy || '', x.occupancy.toFixed(4), x.nEvents, x.meanResidencePs.toFixed(1), x.maxResidencePs.toFixed(1), x.D.toFixed(4), x.bestFrame, x.bestD.toFixed(3), x.bestT1.toFixed(1), gs.G.toFixed(4)]);
    });
  }
}
fs.writeFileSync(path.join(OUT, 'hiyd_t2_v3.csv'), rows.map(r => r.map(v => /[,"]/.test(String(v)) ? '"' + v + '"' : v).join(',')).join('\n') + '\n');
console.log(rows.length - 1, 'contact rows');
