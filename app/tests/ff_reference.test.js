// Compares BondX force-field charges with PDB2PQR reference PQR files (independent implementation).
// For every residue present in both, the charge of each heavy atom plus its hydrogens ("folded" charge)
// must agree, and the residue net charges must agree. Histidine tautomers and other protonation
// states are taken from the reference file (PDB2PQR optimises them), so the comparison tests the
// charge tables, the template matching, termini detection and hydrogen folding.
const fs = require('fs'); const path = require('path');
const B = require('../src/engine.js');
const R = p => fs.readFileSync(path.join(__dirname, 'ref', p), 'utf8');
let fails = 0; const out = [];
function compare(pdb, pqr, ff, label) {
  const S = B.parsePDB(R(pdb));
  // some reference files contain the model twice: keep the first copy (until the atom serial restarts)
  let seen = false; const ql = []; for (const l of R(pqr).split(/\r?\n/)) { if (/^ATOM\s+1\s/.test(l)) { if (seen) break; seen = true; } ql.push(l); }
  // whitespace PQR files may omit chain IDs: give each run of increasing residue numbers the next chain ID of the PDB file
  const chains = [...new Set(S.atoms.filter(a => !a.het).map(a => a.chain))]; let ci = 0, last = -1e9;
  const ql2 = ql.map(l => { if (!/^(ATOM|HETATM)/.test(l)) return l; const t = l.trim().split(/\s+/); if (t.length !== 10) return l; const rs = parseInt(t[4]);
    if (rs < last && !/HOH|WAT/.test(t[3])) ci++; last = rs; t.splice(4, 0, chains[Math.min(ci, chains.length - 1)] || 'A'); return t[0].padEnd(6) + ' ' + t.slice(1).join(' '); });
  const Q = B.parsePDB(ql2.join('\n'), { pqr: true });
  // protonation choices from the reference file
  const overrides = {}; const qres = new Map();
  for (const r of Q.residues) {
    qres.set(r.chain + r.resSeq + (r.iCode || ''), r);
    const rn = r.resName.length === 4 && /^[NC]/.test(r.resName) ? r.resName.slice(1) : r.resName;
    const has = n => r.byName[n] !== undefined;
    let v = null;
    if (['HIS', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP'].includes(rn)) v = has('HD1') && has('HE2') ? 'HIP' : has('HD1') ? 'HID' : 'HIE';
    else if (rn === 'ASP' || rn === 'ASH') v = has('HD2') || rn === 'ASH' ? 'ASH' : 'ASP';
    else if (rn === 'GLU' || rn === 'GLH') v = has('HE2') || rn === 'GLH' ? 'GLH' : 'GLU';
    else if (rn === 'LYS' || rn === 'LYN') v = has('HZ3') ? 'LYS' : 'LYN';
    else if (rn === 'CYX' || rn === 'CYM') v = rn;
    if (v) { const sr = S.residues.find(x => x.chain === r.chain && x.resSeq === r.resSeq && (x.iCode || '') === (r.iCode || '') && x.cls === 'protein'); if (sr) overrides[B.resId(sr)] = v; }
  }
  const P = B.mergeParams(B.DEFAULTS, { ff: { name: ff, overrides, chargeSource: 'ff', foldH: true, missing: 'ignore' } });
  const A = B.assignForceField(S, P);
  // fold reference hydrogens onto heavy atoms
  const qf = new Map();
  for (const a of Q.atoms) { if (a.el === 'H') continue; let q = a.q; for (const j of Q.bonds[a.i]) if (Q.atoms[j].el === 'H') q += Q.atoms[j].q; qf.set(a.chain + a.resSeq + (a.iCode || '') + '|' + a.name, q); }
  let nAtoms = 0, nRes = 0, maxDev = 0; const bad = []; const resDiff = [];
  for (const r of S.residues) {
    if (!(r.cls === 'protein' || r.cls === 'water')) continue;
    const ref = qres.get(r.chain + r.resSeq + (r.iCode || '')); if (!ref) continue;
    const info = A.residues[r.idx]; if (info.missing.length) continue; // incomplete residue in the input
    const heavyRef = ref.atoms.filter(i => Q.atoms[i].el !== 'H').map(i => Q.atoms[i].name);
    const heavyOur = r.atoms.filter(i => S.atoms[i].el !== 'H').map(i => S.atoms[i].name.replace(/^OT2$/, 'OXT').replace(/^OT1$/, 'O'));
    if (heavyRef.length !== heavyOur.length) continue; // reference rebuilt atoms
    nRes++;
    let qnet = 0, qref = 0;
    for (const i of r.atoms) { const a = S.atoms[i]; if (a.el === 'H') continue; const key = r.chain + r.resSeq + (r.iCode || '') + '|' + a.name; const rq = qf.get(key); let my = A.q[i]; for (const j of S.bonds[i]) if (S.atoms[j].el === 'H') my += A.q[j]; qnet += my; if (rq === undefined) continue; qref += rq; nAtoms++;
      const d = Math.abs(rq - my); maxDev = Math.max(maxDev, d); if (d > 2e-4) bad.push(`${r.resName}${r.resSeq}@${a.name}: BondX ${my.toFixed(4)} ref ${rq.toFixed(4)} (${info.tpl})`); }
    if (Math.abs(qnet - qref) > 2e-3) resDiff.push(`${r.resName}${r.resSeq} (${info.tpl}): ${qnet.toFixed(3)} vs ${qref.toFixed(3)}`);
  }
  const ok = bad.length === 0 && resDiff.length === 0 && nAtoms >= 50;
  out.push(`${ok ? 'PASS' : 'FAIL'}  ${label}: ${nRes} residues, ${nAtoms} heavy atoms compared, max |Δq| = ${maxDev.toExponential(2)} e` + (ok ? '' : `\n      ${bad.slice(0, 8).join('\n      ')}\n      residue nets: ${resDiff.slice(0, 8).join('; ')}`));
  return { ok, bad, resDiff };
}
const cases = [
  ['1AFS.pdb', '1AFS_whitespace_ff=AMBER.pqr', 'ff14SB', '1AFS  AMBER'],
  ['1AJJ.pdb', '1AJJ_pH7_AMBER.pqr', 'ff14SB', '1AJJ  AMBER pH 7'],
  ['1BX8.pdb', '1BX8_pH7_AMBER.pqr', 'ff14SB', '1BX8  AMBER pH 7'],
  ['1A1P.pdb', '1A1P_assign-only_whitespace_ff=AMBER.pqr', 'ff14SB', '1A1P  AMBER'],
];
for (const c of cases) { const r = compare(...c); if (!r.ok) fails++; }
// CHARMM: PDB2PQR ships CHARMM22 charges; CHARMM36 differs only for Trp. Count residues that differ.
for (const [pdb, pqr, lab] of [['1AFS.pdb', '1AFS_whitespace_ff=CHARMM.pqr', '1AFS'], ['1AJJ.pdb', '1AJJ_pH7_CHARMM.pqr', '1AJJ'], ['1BX8.pdb', '1BX8_pH7_CHARMM.pqr', '1BX8']]) {
  const r = compare(pdb, pqr, 'charmm36', lab + '  CHARMM36 vs CHARMM22 reference');
  const nonTrp = r.bad.filter(s => !/^TRP/.test(s));
  out[out.length - 1] = (nonTrp.length ? 'FAIL' : 'PASS') + out[out.length - 1].slice(4).split('\n')[0] + `  (differences outside Trp: ${nonTrp.length}; Trp atoms differing as expected for C36: ${r.bad.length - nonTrp.length})` + (nonTrp.length ? '\n      ' + nonTrp.slice(0, 10).join('\n      ') : '');
  if (nonTrp.length) fails++;
}
console.log(out.join('\n'));
if (fails) { console.log(`\n${fails} reference comparison(s) failed`); process.exit(1); } else console.log('\nAll reference comparisons passed.');
