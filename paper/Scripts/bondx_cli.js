#!/usr/bin/env node
/* bondx_cli.js - command-line front end to the BondX Studio engine (engine.js + ff_data.js).
   The Python scripts call this file. Every command prints JSON to stdout.
   Commands
     detect-pdb  <pdb> [--P json]                       T1 detection on a PDB/PQR file
     detect-xyz  <xyz> [--P json]                       T1 detection on an XYZ cluster (materials path)
     assign      <pdb> [--P json] [--mol2 file]          force-field protonation and charges (summary)
     ligand-check <pdb>                                  ligand H count and charge from the heavy-atom structure vs. the PDB FORMUL record
     build-pdb   <pdb> --out dir --name tag [--P json] [--ligand RES:CHAIN:SEQ] [--mol2 file] [--edaBasis b]
                                                        T3a QM model (cluster or QM/MM from P.qm.env); Q-Chem NBO + EDA inputs
     build-xyz   <xyz> --out dir --name tag --X n --A n (--Y n | --ring a,b,c) [--P json] [--edaMethod m --edaBasis b]
                                                        two-fragment model for a gas-phase dimer; Q-Chem NBO + EDA inputs
     gate        <manifest.json> <nbo_output> [--threshold t] [--ry 0|1]
     eda         <eda_output>
     bxqi        --json '{"G":..,"E2":..,"eda":{..}}'
   Indices for build-xyz are 1-based positions in the XYZ file. */
'use strict';
const fs = require('fs'), path = require('path');
const B = require('./engine.js');
const argv = process.argv.slice(2); const cmd = argv[0];
const opt = (k, d) => { const i = argv.indexOf('--' + k); return i >= 0 ? argv[i + 1] : d; };
const P = B.mergeParams(B.defaults(), JSON.parse(opt('P', '{}')));
const out = o => process.stdout.write(JSON.stringify(o, null, 1) + '\n');
const read = f => fs.readFileSync(f, 'utf8');
const loadPDB = f => B.parsePDB(read(f), { pqr: /\.pqr$/i.test(f), altloc: P.struct.altloc, bondTol: P.struct.bondTol });
const slimHit = (S, h) => ({ key: h.key, donor: h.donorLabel, donorResName: S.residues[h.donorRes].resName, donorRes: h.donorRes, X: h.Xname, A: h.Aname, halogen: h.halogen,
  acceptor: h.accLabel, accResName: S.residues[h.accRes].resName, accResSeq: S.residues[h.accRes].resSeq, accChain: S.residues[h.accRes].chain, accAtom: h.accName, accEl: h.accEl, accClass: h.accClass, type: h.type,
  lewis: h.lewis, d: h.d, sumVdw: h.sumVdw, R: h.R, theta1: h.theta1, theta2: h.theta2, gd: h.gd, gt: h.gt, g2: h.g2, G: h.G, status: h.status });

// Level-of-theory names per program. [method keyword, dispersion: none | bj | zero]
// ORCA has no wB97X-D (Chai/Head-Gordon); its nearest analogue wB97X-D3 is used and labelled as such.
const METH = {
  'wB97X-D':    { g: ['wB97XD', 'none'],  o: ['wB97X-D3', 'none'], q: ['wB97X-D', 'none'] },
  'M06-2X':     { g: ['M062X', 'none'],   o: ['M062X', 'none'],    q: ['M06-2X', 'none'] },
  'B3LYP-D3BJ': { g: ['B3LYP', 'bj'],     o: ['B3LYP', 'bj'],      q: ['B3LYP', 'bj'] },
  'PBE0-D3BJ':  { g: ['PBE1PBE', 'bj'],   o: ['PBE0', 'bj'],       q: ['PBE0', 'bj'] },
  'HF':         { g: ['HF', 'none'],      o: ['HF', 'none'],       q: ['HF', 'none'] },
  'wB97X-V':    { g: null,                o: ['wB97X-V', 'none'],  q: ['wB97X-V', 'none'] }
};
// Gaussian 16 has no built-in def2-TZVPD/def2-SVPD; those jobs are written for ORCA and Q-Chem only.
const BAS = { 'def2-SVP': { g: 'def2SVP', o: 'def2-SVP', q: 'def2-SVP' }, 'def2-TZVP': { g: 'def2TZVP', o: 'def2-TZVP', q: 'def2-TZVP' },
  'def2-TZVPD': { g: null, o: 'def2-TZVPD', q: 'def2-TZVPD' }, 'def2-SVPD': { g: null, o: 'def2-SVPD', q: 'def2-SVPD' }, 'def2-QZVP': { g: 'def2QZVP', o: 'def2-QZVP', q: 'def2-QZVP' } };
const PROGS = (opt('programs', 'gaussian,orca,qchem')).split(',');
function writeInputs(model, outDir, name, pcs, opts) {
  const files = []; const w = (sub, f, txt) => { const d = path.join(outDir, sub); fs.mkdirSync(d, { recursive: true }); fs.writeFileSync(path.join(d, f), txt); files.push(sub + '/' + f); };
  const M = METH[opts.method], Bs = BAS[opts.basis]; if (!M || !Bs) throw new Error('unknown level ' + opts.method + '/' + opts.basis);
  const memGB = Math.round(opts.memMB / 1000), np = opts.nproc;
  if (PROGS.includes('gaussian') && M.g && Bs.g) w('gaussian', `${name}.gjf`, B.writeGaussian(model, { job: name, nproc: np, mem: memGB + 'GB', method: M.g[0], disp: M.g[1], basis: Bs.g, solvent: opts.solvent, eps: opts.eps, nbo7: false }, pcs));
  if (PROGS.includes('orca') && M.o && Bs.o) { const r = B.writeORCA(model, { job: name, nproc: np, maxcore: Math.floor(opts.memMB * 0.75 / np), method: M.o[0], disp: M.o[1], basis: Bs.o, solvent: opts.solvent, eps: opts.eps }, pcs);
    w('orca', `${name}.inp`, r.inp); if (r.pc) w('orca', `${name}.pc`, r.pc); }
  if (PROGS.includes('qchem') && M.q && Bs.q) w('qchem', `${name}.in`, B.writeQChemNBO(model, { method: M.q[0], disp: M.q[1], basis: Bs.q, memMB: opts.memMB, solvent: opts.solvent, eps: opts.eps }, pcs));
  if (opts.eda !== false && !pcs) { const E = METH[opts.edaMethod], EB = BAS[opts.edaBasis];
    w('qchem_eda', `${name}_eda.in`, B.writeQChemEDA(model, { method: E.q[0], disp: E.q[1], basis: EB.q, eda2: 1, bsse: !opts.noBSSE, memMB: opts.memMB, solvent: opts.edaSolvent, eps: opts.eps })); }
  // Interaction-energy jobs for Gaussian and ORCA, run next to every Q-Chem ALMO-EDA job.
  // Gaussian and ORCA have no ALMO-EDA; these give the counterpoise-corrected interaction energy at a matching level,
  // which is compared with the ALMO-EDA total.
  if (opts.eda !== false && !pcs) {
    const m = model.manifest, fq = m.fragCharge, Zs = B.Z;
    const ne = f => model.qm.filter(a => (a.frag === 1) === (f === 1)).reduce((t, a) => t + (Zs[a.el] || 0), 0) - (fq[f] || 0);
    const mult = f => (ne(f) % 2 ? 2 : 1); const xyz = (a, sym) => `${sym.padEnd(14)} ${a.x.toFixed(6).padStart(12)} ${a.y.toFixed(6).padStart(12)} ${a.z.toFixed(6).padStart(12)}`;
    const el = a => a.el[0] + a.el.slice(1).toLowerCase();
    const memGB = Math.round(opts.memMB / 1000);
    if (PROGS.includes('gaussian')) {
      const G = ['%nprocshared=' + opts.nproc, '%mem=' + memGB + 'GB', `%chk=${name}_cp.chk`, `#p wB97XD/def2TZVP counterpoise=2 nosymm scf=(tight,xqc) int=ultrafine`, '',
        `BondX ${m.id} | counterpoise interaction energy (pairs with Q-Chem ALMO-EDA)`, '', `${m.totalCharge} ${m.multiplicity} ${fq[1] || 0} ${mult(1)} ${fq[2] || 0} ${mult(2)}`];
      model.qm.forEach(a => G.push(xyz(a, `${el(a)}(Fragment=${a.frag === 1 ? 1 : 2})`))); G.push('', '');
      w('gaussian_cp', `${name}_cp.gjf`, G.join('\n'));
    }
    if (PROGS.includes('orca')) {
      const job = (tag, q, mu, ghostOf) => { const L = [`# BondX ${m.id} - counterpoise job ${tag} (E_int = E_AB - E_A(AB basis) - E_B(AB basis))`,
        `! wB97X-V def2-TZVPD RIJCOSX def2/J TightSCF${opts.cpRobust ? ' DefGrid3 SlowConv' : ''}`, `%pal nprocs ${opts.nproc} end`, `%maxcore ${Math.floor(opts.memMB * 0.75 / opts.nproc)}`,
        ...(opts.cpRobust ? ['%scf', '  Sthresh 1e-6    # remove near-linearly-dependent diffuse functions (same basis in AB, A and B)', '  MaxIter 500', 'end'] : []), '', `* xyz ${q} ${mu}`];
        model.qm.forEach(a => { const f = a.frag === 1 ? 1 : 2; L.push(xyz(a, el(a) + (ghostOf && ghostOf !== f ? ':' : ''))); }); L.push('*', ''); return L.join('\n'); };
      w('orca_cp', `${name}_cp_AB.inp`, job('AB', m.totalCharge, m.multiplicity, null));
      w('orca_cp', `${name}_cp_A.inp`, job('A (fragment 1, fragment 2 as ghost atoms)', fq[1] || 0, mult(1), 1));
      w('orca_cp', `${name}_cp_B.inp`, job('B (fragment 2, fragment 1 as ghost atoms)', fq[2] || 0, mult(2), 2));
    }
  }
  w('manifests', `${name}_manifest.json`, JSON.stringify(model.manifest, null, 1));
  w('xyz', `${name}.xyz`, B.writeXYZ(model));
  if (pcs) w('mm_charges', `${name}_mm_charges.csv`, B.chargesCSV(pcs));
  return files;
}
const QOPT = () => ({ method: opt('method', 'wB97X-D'), basis: opt('basis', 'def2-TZVP'), memMB: +opt('mem', 32000), nproc: +opt('nproc', 16),
  solvent: P.qm.env === 'cluster' && P.qm.pcm, eps: P.qm.eps, edaMethod: opt('edaMethod', 'wB97X-V'), edaBasis: opt('edaBasis', 'def2-TZVPD'), edaSolvent: false, eda: opt('noEda') ? false : true, noBSSE: !!opt('noBSSE'), cpRobust: !!opt('cpRobust') });

if (cmd === 'detect-pdb') {
  const t0 = process.hrtime.bigint(); const S = loadPDB(argv[1]); const t1 = process.hrtime.bigint();
  const { hits } = B.detectXB(S, P); const t2 = process.hrtime.bigint();
  out({ file: path.basename(argv[1]), atoms: S.atoms.length, residues: S.residues.length, hasH: S.hasH,
    halogenatedLigands: S.residues.filter(r => r.cls === 'ligand' && r.halogens.some(x => x !== 'F')).map(r => ({ resName: r.resName, chain: r.chain, resSeq: r.resSeq, halogens: r.halogens })),
    msParse: Number(t1 - t0) / 1e6, msDetect: Number(t2 - t1) / 1e6, hits: hits.map(h => slimHit(S, h)) });
} else if (cmd === 'detect-xyz') {
  const M = B.parseXYZ(read(argv[1])); out({ atoms: M.atoms.length, hits: B.detectMaterialXB(M, P).map(h => ({ X: h.X + 1, A: h.A + 1, Y: h.Y + 1, Xlabel: h.Xlabel, Ylabel: h.Ylabel, halogen: h.halogen, accEl: h.accEl, d: h.d, theta1: h.theta1, G: h.G, gd: h.gd, gt: h.gt, status: h.status })) });
} else if (cmd === 'ligand-check') {
  const S = loadPDB(argv[1]); const rows = [];
  for (const R of S.residues) {
    if (R.cls !== 'ligand') continue;
    const lp = B.ligandProtonation(S, S.frames[0], R, 'rules'); const fileH = R.atoms.filter(i => S.atoms[i].el === 'H').length;
    const fc = lp.metal ? null : B.formulCheck(S, R, lp.nHtotal, lp.charge);
    rows.push({ resName: R.resName, chain: R.chain, resSeq: R.resSeq, heavy: R.atoms.length - fileH, fileH, metal: !!lp.metal, modelH: lp.nHtotal, charge: lp.charge,
      formul: !!fc, complete: fc ? fc.complete : null, links: fc ? fc.links : null, expH: fc ? fc.expH : null, hOK: fc ? fc.hOK : null, notes: lp.notes.join('; ') });
  }
  out({ file: path.basename(argv[1]), ligands: rows });
} else if (cmd === 'assign') {
  const S = loadPDB(argv[1]); const mol2 = opt('mol2') ? B.mol2ForStructure(B.parseMOL2(read(opt('mol2')), true), S) : null;
  const ffa = B.assignForceField(S, P, { mol2 });
  out({ ff: ffa.ff, water: ffa.water, total: ffa.total, usedPQR: ffa.usedPQR, warnings: ffa.warnings,
    residues: ffa.residues.map(r => ({ id: r.id, resName: r.resName, cls: r.cls, tpl: r.tpl, why: r.why, net: r.net, formal: r.formal, missing: r.missing, extra: r.extra, foldedH: r.foldedH, source: r.source })),
    titratable: B.titratableResidues(S, P, ffa).map(t => ({ id: t.id, resName: t.resName, variant: t.variant, choices: t.choices, net: t.net, why: t.why })) });
} else if (cmd === 'build-pdb') {
  const t0 = Date.now(); const S = loadPDB(argv[1]); const mol2 = opt('mol2') ? B.mol2ForStructure(B.parseMOL2(read(opt('mol2')), true), S) : null;
  const ffa = B.assignForceField(S, P, { mol2 }); const { hits } = B.detectXB(S, P);
  let xb = hits.filter(h => h.status === 'XB');
  const lig = opt('ligand'); if (lig) { const [rn, ch, sq] = lig.split(':'); xb = xb.filter(h => { const R = S.residues[h.donorRes]; return R.resName === rn && (!ch || R.chain === ch) && (!sq || String(R.resSeq) === sq); }); }
  if (!xb.length) { out({ error: 'no XB for ligand ' + (lig || '(any)') }); process.exit(0); }
  const donor = xb[0].donorRes; xb = xb.filter(h => h.donorRes === donor).map((h, k) => Object.assign({ cid: 'XB' + (k + 1) }, h));
  const model = B.buildQMModel(S, 0, xb, { P, ffa, mol2, sourceName: path.basename(argv[1]) });
  let pcs = null; if (P.qm.env === 'qmmm') pcs = B.embeddingCharges(S, 0, model.include, { P, ffa, qmResidueCharges: model.qmResCharge });
  const mf = model.manifest; if (pcs) mf.embedding = { source: pcs.source, water: pcs.water, n: pcs.charges.length, net: pcs.net, radius: pcs.radius, boundary: pcs.boundary, M1: pcs.M1, M2: pcs.M2, cut: pcs.cut, notes: pcs.notes };
  const o = QOPT(); mf.level = `${o.method}/${o.basis}`; mf.parameters = { qm: P.qm, mm: P.mm, ff: P.ff };
  const files = writeInputs(model, opt('out'), opt('name'), pcs, o);
  const e = model.qm.reduce((s, a) => s + (B.Z[a.el] || 0), 0) - mf.totalCharge;
  out({ name: opt('name'), id: mf.id, nAtoms: mf.nAtoms, addedH: mf.addedH, linkH: mf.linkH, fragCharge: mf.fragCharge, totalCharge: mf.totalCharge, multiplicity: mf.multiplicity, electrons: e,
    nHeavy: model.qm.filter(a => a.el !== 'H').length, xbs: mf.xbs.map(x => ({ id: x.id, label: x.label, d: x.d, theta1: x.theta1, halogen: x.halogen })), residueCharges: mf.residueCharges, warnings: mf.warnings,
    embedding: mf.embedding || null, files, ms: Date.now() - t0 });
} else if (cmd === 'build-xyz') {
  const M = B.parseXYZ(read(argv[1])); const X = +opt('X') - 1, A = +opt('A') - 1; const ring = opt('ring') ? opt('ring').split(',').map(n => +n - 1) : null; const Y = ring ? ring[0] : +opt('Y') - 1;
  const at = M.atoms; const d = B.dist(at[X].c, at[Y].c);
  let Yc = at[Y].c; if (ring) { Yc = [0, 1, 2].map(k => ring.reduce((s, i) => s + at[i].c[k], 0) / ring.length); }
  const hit = { cid: 'XB1', X, A, Y, Ypos: at[Y].c, image: '0 0 0', halogen: at[X].el, Xlabel: at[X].label, Ylabel: ring ? 'ring' : at[Y].label, d: B.dist(at[X].c, Yc), theta1: B.angle(at[A].c, at[X].c, Yc) };
  const model = B.buildMaterialQM(M, hit, { P, sourceName: path.basename(argv[1]) });
  const mf = model.manifest; const map = new Map(); mf.atoms.forEach((a, k) => map.set(a.src.replace(/'$/, ''), k + 1));
  if (ring) { const xb = mf.xbs[0]; xb.type = 'pi'; xb.ring = ring.map(i => map.get(at[i].label)); xb.Y = null; xb.label = `${at[X].label} ··· ring(${ring.map(i => at[i].label).join(',')})`; }
  mf.xbs[0].d = hit.d; mf.xbs[0].theta1 = hit.theta1;
  const o = QOPT(); o.solvent = false; mf.level = `${o.method}/${o.basis}`;
  const files = writeInputs(model, opt('out'), opt('name'), null, o);
  out({ name: opt('name'), id: mf.id, nAtoms: mf.nAtoms, frag1: mf.atoms.filter(a => a.frag === 1).length, frag2: mf.atoms.filter(a => a.frag === 2).length, xb: mf.xbs[0], warnings: mf.warnings, files });
} else if (cmd === 'gate') {
  const mf = JSON.parse(read(argv[1])); const nbo = B.parseNBO(read(argv[2]));
  if (nbo.error) { out({ error: nbo.error }); process.exit(0); }
  const ver = B.verifyGeometry(mf, nbo.geometry, P.nbo.geomTol);
  const g = B.nboGate(mf, nbo, { threshold: opt('threshold') !== undefined ? +opt('threshold') : P.nbo.threshold, includeRY: opt('ry') !== undefined ? opt('ry') === '1' : P.nbo.includeRY });
  out({ program: nbo.program, terminated: nbo.terminated, nEntries: nbo.entries.length, printThreshold: nbo.threshold, verification: ver, fragNPA: B.fragmentNPA(mf, nbo),
    gate: g.map(x => ({ id: x.id, label: x.label, ionicX: x.ionicX, E2_core: x.E2_core, E2_total: x.E2_total, contrib: x.contrib, verified: x.verified, npaX: x.npaX, npaY: x.npaY, nRows: x.rows.length })) });
} else if (cmd === 'eda') {
  const r = B.parseEDA(read(argv[1])); out(Object.assign(r, { descriptors: r.error ? null : B.edaDescriptors(r, P) }));
} else if (cmd === 'gscore') {
  // geometric score for one contact: --el I --acc O|N|S|AR --d 3.1 --t1 179 [--t2 120] [--alpha 10]
  const el = opt('el').toUpperCase(), acc = opt('acc').toUpperCase(), d = +opt('d'), t1 = +opt('t1');
  const isPi = acc === 'AR' || acc === 'PI'; const c = (P.geom.cutoffs[el] || {})[isPi ? 'AR' : acc];
  const g = c ? B.geometricScore(el, isPi ? 'AR' : acc, d, t1, opt('t2') !== undefined ? +opt('t2') : null, P, isPi, opt('alpha') !== undefined ? +opt('alpha') : null) : null;
  const pass = !!c && d <= c && t1 >= P.geom.theta1Min && t1 <= P.geom.theta1Max && (!isPi || opt('alpha') === undefined || +opt('alpha') <= P.geom.piNormalMax);
  out({ covered: !!c, cutoff: c || null, pass, G: g && pass ? g.G : 0, gd: g ? g.gd : null, gt: g ? g.gt : null, g2: g ? g.g2 : null, R: g ? g.R : null });
} else if (cmd === 'bxqi') {
  out(B.computeBXQI(JSON.parse(opt('json')), P));
} else {
  process.stderr.write(fs.readFileSync(__filename, 'utf8').split('\n').slice(1, 17).join('\n') + '\n'); process.exit(1);
}
