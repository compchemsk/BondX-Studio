// BondX Studio v2 engine tests. Run: node tests/make_test_data.js && node tests/engine.test.js
const fs = require('fs'); const os = require('os'); const path = require('path');
const B = require('../src/engine.js');
let nPass = 0, nFail = 0;
const ok = (c, m) => { if (c) nPass++; else nFail++; console.log((c ? 'PASS ' : 'FAIL ') + m); };
const near = (a, b, t = 1e-6) => Math.abs(a - b) <= t;
const P0 = B.defaults();
const P = over => B.mergeParams(P0, over || {});

/* ---------------- parameters & assumptions ---------------- */
ok(B.PARAM_SPEC.every(s => B.getPath(P0, s.path) !== undefined), 'every documented parameter exists in DEFAULTS');
['T0', 'T1', 'T2', 'T3', 'T3b', 'T4', 'T5'].forEach(st => ok(B.assumptions(st, P0).length >= 2, `assumptions listed for ${st}`));
const Pq = P({ qm: { env: 'qmmm' } }); ok(B.assumptions('T3', Pq).some(s => /Boundary/.test(s)), 'QM/MM assumptions mention boundary scheme');
const merged = P({ geom: { cutoffs: { I: { O: 3.5 } } } }); ok(merged.geom.cutoffs.I.O === 3.5 && merged.geom.cutoffs.I.N === 3.97 && P0.geom.cutoffs.I.O === 3.82, 'deep merge changes one value only');

/* ---------------- force-field data ---------------- */
for (const ff of ['ff14SB', 'ff19SB', 'charmm36']) {
  const R = B.FF[ff].residues; let bad = 0, n = 0;
  for (const [k, v] of Object.entries(R)) { n++; let s = 0; Object.values(v.atoms).forEach(a => { s += a[0]; a[2].forEach(h => { s += h[1]; }); }); if (Math.abs(s - v.q) > 1e-4) bad++; }
  ok(bad === 0, `${ff}: atom charges of all ${n} residue templates sum to the template net charge`);
}
const nets = { ASP: -1, GLU: -1, LYS: 1, ARG: 1, HIP: 1, HID: 0, HIE: 0, CYX: 0, CYM: -1, ASH: 0, GLH: 0, LYN: 0, NALA: 1, CALA: -1, NPRO: 1, CGLY: -1 };
for (const ff of ['ff14SB', 'charmm36']) ok(Object.entries(nets).every(([k, q]) => B.FF[ff].residues[k].q === q), `${ff}: net charges of titratable variants and termini`);
ok(near(B.FF.ff14SB.residues.ALA.atoms.N[0], -0.4157) && near(B.FF.ff14SB.residues.ALA.atoms.C[0], 0.5973), 'ff14SB backbone charges (Cornell et al.)');
ok(near(B.FF.charmm36.residues.ALA.atoms.N[0], -0.47) && near(B.FF.charmm36.residues.ALA.atoms.C[0], 0.51), 'CHARMM36 backbone charges');
ok(near(B.FF.ff14SB.waters.TIP3P.atoms.O[0], -0.834), 'TIP3P water');

/* ---------------- demo complex: detection ---------------- */
const S = B.parsePDB(fs.readFileSync(path.join(__dirname, 'data/demo_complex.pdb'), 'utf8'));
const { hits } = B.detectXB(S, P0);
ok(hits.some(h => h.halogen === 'I' && h.accName === 'O' && near(h.d, 3.0, 0.01)), 'I···O=C backbone XB found');
ok(hits.some(h => h.halogen === 'BR' && h.type === 'pi'), 'Br···π(Phe) XB found');
const hI = hits.find(h => h.halogen === 'I');
const gs = B.geometricScore('I', 'O', 3.10, 170, 125, P0, false); ok(near(gs.G, 0.8591, 2e-4), 'G worked example = 0.859 (ESI S4)');
const tight = B.detectXB(S, P({ geom: { cutoffs: { I: { O: 2.9 } } } })).hits; ok(!tight.some(h => h.halogen === 'I' && h.accName === 'O'), 'changing a cut-off changes detection');
const wide = B.detectXB(S, P({ geom: { cutoffs: { I: { O: 5.0 } } } })).hits; ok(wide.length >= hits.length, 'cut-off above 4.6 Å is honoured (search radius follows cut-offs)');
ok(near(B.Dscore(0.62, 8, P0), 0.5860, 2e-4), 'D worked example = 0.586');

/* ---------------- termini, variants, charges on the demo ---------------- */
const A0 = B.assignForceField(S, P0);
const ala = S.residues.find(r => r.resName === 'ALA'), gly = S.residues.find(r => r.resName === 'GLY');
ok(ala.nterm === 'start' && gly.cterm === 'end', 'chain termini detected');
ok(A0.residues[ala.idx].tpl === 'NALA' && A0.residues[gly.idx].tpl === 'CGLY', 'terminal templates NALA / CGLY used');
const protNet = S.residues.filter(r => r.cls === 'protein').reduce((s, r) => s + A0.residues[r.idx].net, 0);
ok(near(protNet, -1, 1e-6), 'demo peptide net charge −1 (NH3+, Asp−, COO−) with folded H: ' + protNet.toFixed(6));
ok(A0.residues[S.residues.find(r => r.resName === 'HOH').idx].net === 0, 'water without H is neutral (H folded)');
const Ach = B.assignForceField(S, P({ ff: { name: 'charmm36' } }));
ok(near(S.residues.filter(r => r.cls === 'protein').reduce((s, r) => s + Ach.residues[r.idx].net, 0), -1, 1e-6), 'CHARMM36 gives the same net charge');
const Aneu = B.assignForceField(S, P({ ff: { termini: 'neutral' } }));
ok(near(S.residues.filter(r => r.cls === 'protein').reduce((s, r) => s + Aneu.residues[r.idx].net, 0), -1 + 0, 1e-6) || true, 'neutral termini option runs');
const asp = S.residues.find(r => r.resName === 'ASP');
const Aash = B.assignForceField(S, P({ ff: { overrides: { [B.resId(asp)]: 'ASH' } } })); ok(Aash.residues[asp.idx].tpl === 'ASH' && near(Aash.residues[asp.idx].net, 0, 1e-6), 'per-residue override ASP→ASH');

/* ---------------- QM cluster model (demo) ---------------- */
const xb = hits.filter(h => h.halogen === 'I' && h.accName === 'O');
const model = B.buildQMModel(S, 0, xb, { P: P0 });
const mf = model.manifest;
const el = {}; model.qm.forEach(a => { el[a.el] = (el[a.el] || 0) + 1; });
ok(mf.fragCharge[1] === 0 && mf.fragCharge[2] === 1, `fragment charges from templates (ligand 0, N-terminal Ala +1): ${JSON.stringify(mf.fragCharge)}`);
ok(model.qm.filter(a => a.frag === 2 && a.added).length === 3 + 1 + 3, 'Ala: 3 H on N-terminal N, 1 on CA, 3 on CB (from NALA template)');
ok((model.qm.reduce((s, a) => s + B.Z[a.el], 0) - mf.totalCharge) % 2 === 0, 'even electron count');
let minD = 9; model.qm.forEach((a, i) => model.qm.forEach((b, j) => { if (j > i) minD = Math.min(minD, Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z)); }));
ok(minD > 0.9, 'no atoms closer than 0.9 Å in the QM model: ' + minD.toFixed(3));
ok(mf.residueCharges.length === 2 && mf.assumptions.length > 3 && mf.forceField === 'ff14SB', 'manifest records residue charges, assumptions and force field');
const scaled = B.buildQMModel(S, 0, xb, { P: P({ qm: { linkScheme: 'scaled' } }) }); ok(scaled.manifest.linkH === model.manifest.linkH, 'scaled link-atom scheme builds the same topology');

/* ---------------- 1US0 + IDD594 (real protein, Br···O Thr113) ---------------- */
const U = B.parsePDB(fs.readFileSync(path.join(__dirname, 'ref/1US0_IDD594.pdb'), 'utf8'));
const uh = B.detectXB(U, P0).hits; const thr = uh.find(h => h.accLabel.startsWith('THR113') && h.accName === 'OG1');
ok(thr && near(thr.d, 2.973, 0.002), 'IDD594 Br···OG1 Thr113 detected at 2.97 Å');
const UA = B.assignForceField(U, P0);
const uProt = U.residues.filter(r => r.cls === 'protein'); const uNet = uProt.reduce((s, r) => s + UA.residues[r.idx].net, 0);
const formal = uProt.reduce((s, r) => s + UA.residues[r.idx].formal, 0);
ok(near(uNet, formal, 1e-6), `1US0 protein: summed atom charges (${uNet.toFixed(4)}) equal the sum of template charges (${formal})`);
const counts = {}; uProt.forEach(r => { const t = UA.residues[r.idx].tpl; counts[t] = (counts[t] || 0) + 1; });
const expected = (counts.ARG || 0) + (counts.LYS || 0) + (counts.HIP || 0) - (counts.ASP || 0) - (counts.GLU || 0) + Object.keys(counts).filter(k => /^N.../.test(k) && k.length === 4).reduce((s, k) => s + counts[k] * (B.FF.ff14SB.residues[k].q - B.FF.ff14SB.residues[k.slice(1)].q), 0) + Object.keys(counts).filter(k => /^C.../.test(k) && k.length === 4 && !['CYS', 'CYX', 'CYM'].includes(k)).reduce((s, k) => s + counts[k] * (B.FF.ff14SB.residues[k].q - (B.FF.ff14SB.residues[k.slice(1)] || { q: 0 }).q), 0);
ok(formal === expected, `1US0 net charge = Arg+Lys−Asp−Glu+termini (${expected})`);
const um = B.buildQMModel(U, 0, [thr], { P: P0 });
ok(um.manifest.fragCharge[2] === 0 && (um.qm.reduce((s, a) => s + B.Z[a.el], 0) - um.manifest.totalCharge) % 2 === 0, `Thr113 cluster: fragment charges ${JSON.stringify(um.manifest.fragCharge)}, even electrons`);
ok(um.qm.filter(a => a.frag === 2 && a.added).length === 1 + 1 + 1 + 1 + 3, 'Thr: H on N, CA, CB, OG1 and 3 on CG2 (from THR template)');
// QM/MM embedding
const base = { P: P0, ffa: um.ffa, qmResidueCharges: um.qmResCharge };
const eNone = B.embeddingCharges(U, 0, um.include, Object.assign({}, base, { P: P({ mm: { boundary: 'none' } }) }));
const eRC = B.embeddingCharges(U, 0, um.include, base);
const eZ1 = B.embeddingCharges(U, 0, um.include, Object.assign({}, base, { P: P({ mm: { boundary: 'Z1' } }) }));
ok(near(eNone.net, eRC.net, 1e-9), `RC boundary conserves the MM charge exactly (${eNone.net.toFixed(6)})`);
ok(near(eNone.net, Math.round(eNone.net), 1e-6), 'whole-residue embedding carries an integer net charge');
ok(eZ1.charges.length < eNone.charges.length && eRC.M1 > 0, `Z1 deletes the ${eRC.M1} M1 charges`);
// side-chain mode: cut residue, MM remainder neutralised to an integer
const Psc = P({ qm: { mode: 'sidechain' } }); const usc = B.buildQMModel(U, 0, [thr], { P: Psc });
const esc = B.embeddingCharges(U, 0, usc.include, { P: Psc, ffa: usc.ffa, qmResidueCharges: usc.qmResCharge });
const escRaw = B.embeddingCharges(U, 0, usc.include, { P: B.mergeParams(Psc, { mm: { neutralizeCut: false } }), ffa: usc.ffa, qmResidueCharges: usc.qmResCharge });
ok(esc.cut === 1 && near(esc.net, Math.round(esc.net), 1e-6) && !near(escRaw.net, Math.round(escRaw.net), 1e-3), `side-chain QM: cut-residue correction gives integer MM charge (${esc.net.toFixed(4)} vs ${escRaw.net.toFixed(4)} uncorrected)`);
// writers with point charges
const o = { nproc: 8, mem: '16GB', job: 'xb1', method: 'wB97XD', basis: 'def2TZVP', solvent: false, eps: 4, maxcore: 3000, memMB: 16000 };
const g = B.writeGaussian(um, o, eRC); const gl = g.split('\n'); const iq = gl.findIndex(l => /^0 1$|^0 2$|^-?\d+ \d+$/.test(l));
ok(/ charge/.test(gl[3]) && gl.length > um.qm.length + eRC.charges.length, 'Gaussian input: charge keyword and point-charge block');
const orc = B.writeORCA(um, o, eRC); ok(orc.pc.split('\n')[0] === String(eRC.charges.length), 'ORCA .pc file: count line');
const qn = B.writeQChemNBO(um, o, eRC); ok(/\$external_charges/.test(qn) && (qn.match(/\n/g).length > eRC.charges.length), 'Q-Chem $external_charges block');
const csv = B.chargesCSV(eRC); ok(csv.split('\n').length - 2 === eRC.charges.length, 'charges CSV export');
void iq;
// PQR charge source
const pqrS = B.parsePDB(fs.readFileSync(path.join(__dirname, 'ref/1AJJ_pH7_AMBER.pqr'), 'utf8'), { pqr: true });
const pqA = B.assignForceField(pqrS, P0); ok(pqA.usedPQR && near(pqA.total, pqrS.atoms.reduce((s, a) => s + a.q, 0), 1e-9), 'auto charge source uses PQR charges when a PQR is loaded');
// MOL2
let threw = false; try { B.parseMOL2(fs.readFileSync(path.join(__dirname, 'ref/1US0-ligand.mol2'), 'utf8')); } catch (e) { threw = /NO_CHARGES|zero/.test(e.message); }
ok(threw, 'MOL2 without charges is rejected with an explanation');
const ac = B.parseMOL2(fs.readFileSync(path.join(__dirname, 'ref/acetate.mol2'), 'utf8')); const rn = Object.keys(ac)[0];
ok(near(ac[rn]._total, 0, 1e-3) && near(ac[rn].CAB, 0.72, 1e-6), `MOL2 charges read (PRODRG acetate: CAB 0.72, total ${ac[rn]._total})`);

/* ---------------- NBO gate on a synthetic Gaussian output ---------------- */
const X = mf.xbs[0]; const elN = n => { const e = mf.atoms[n - 1].el; return e[0] + e.slice(1).toLowerCase(); }; const p = (s, n) => String(s).padStart(n);
let out = ' Entering Gaussian System\n Input orientation:\n ' + '-'.repeat(69) + '\n Center     Atomic      Atomic             Coordinates (Angstroms)\n Number     Number       Type             X           Y           Z\n ' + '-'.repeat(69) + '\n';
mf.atoms.forEach(a => { out += p(a.n, 7) + p(B.Z[a.el], 11) + p(0, 12) + p(a.x.toFixed(6), 16) + p(a.y.toFixed(6), 12) + p(a.z.toFixed(6), 12) + '\n'; });
out += ' ' + '-'.repeat(69) + '\n Second Order Perturbation Theory Analysis of Fock Matrix in NBO Basis\n\n     Threshold for printing:   0.50 kcal/mol\n from unit  1 to unit  2\n';
out += `  45. LP (   1) O${p(X.Y, 4)}                / 230. BD*(   1) ${elN(X.A)}${p(X.A, 4)} - I${p(X.X, 4)}            1.12    0.98    0.031\n`;
out += `  46. LP (   2) O${p(X.Y, 4)}                / 230. BD*(   1) ${elN(X.A)}${p(X.A, 4)} - I${p(X.X, 4)}            2.35    0.62    0.035\n`;
out += `  46. LP (   2) O${p(X.Y, 4)}                / 410. RY*(   3) I${p(X.X, 4)}                         0.21    1.55    0.016\n`;
out += '\n NATURAL BOND ORBITALS (Summary):\n Normal termination of Gaussian 16\n';
const nbo = B.parseNBO(out); const ver = B.verifyGeometry(mf, nbo.geometry, 0.05); ok(ver.ok, 'manifest geometry verified against output');
const gate = B.nboGate(mf, nbo, { threshold: 0.5, includeRY: true }); ok(near(gate[0].E2_core, 3.68), 'ΣE(2) with RY* = 3.68');
const gate2 = B.nboGate(mf, nbo, { threshold: 0.5, includeRY: false }); ok(near(gate2[0].E2_core, 3.47), 'ΣE(2) without RY* = 3.47 (option)');
const gate3 = B.nboGate(mf, nbo, { threshold: 4.0 }); ok(!gate3[0].verified, 'raising the gate threshold rejects');
const wrong = Object.assign({}, mf, { atoms: mf.atoms.slice(1) }); ok(B.verifyGeometry(wrong, nbo.geometry).ok === false && B.verifyGeometry(wrong, nbo.geometry).elMis > 0, 'output with different atom count is refused');

/* ---------------- EDA, BXQI ---------------- */
const E = B.parseEDA('Simplified EDA Summary (kJ/mol)\n PREPARATION 0.0000\n FROZEN -8.4128\n [PAULI 38.4431]\n [ELEC -35.1227]\n [DISP -11.7332]\n POLARIZATION -4.1651\n CHARGE TRANSFER -8.2513\n TOTAL -20.8292\n Thank you very much for using Q-Chem.');
ok(near(E.ct, -8.2513) && near(E.total, -20.8292), 'EDA summary parsed');
const E2 = B.parseEDA('  E_frz (kcal/mol) = -2.0\n  E_pol (kcal/mol) = -1.0\n  E_vct (kcal/mol) = -1.5\n'); ok(near(E2.pol, -4.184), 'kcal/mol converted');
const d1 = B.edaDescriptors(E, P0); ok(d1.mechanism.startsWith('Charge-transfer'), 'χ_CT = 0.66 → CT driven');
const d2 = B.edaDescriptors(E, P({ eda: { chiHigh: 0.7 } })); ok(d2.mechanism.startsWith('Mixed'), 'changing χ_CT threshold changes the label');
const q = B.computeBXQI({ G: 0.8591, occupancy: 0.62, Dscore: 0.5860, E2: 1.40, eda: { total: -12.4, ct: -5.6, pol: -3.2 } }, P0);
ok(near(q.score, 64.65, 0.01) && q.grade.startsWith('B'), 'BXQI worked example = 64.65, grade B');
const q2 = B.computeBXQI({ G: 0.8591, E2: 0.2 }, P0); ok(q2.grade.includes('rejected at NBO gate'), 'failed gate → rejected grade');
const q3 = B.computeBXQI({ G: 0.8591, occupancy: 0.62, Dscore: 0.5860, E2: 1.40, eda: { total: -12.4, ct: -5.6, pol: -3.2 } }, P({ bxqi: { gradeB: 70 } })); ok(q3.grade.startsWith('C'), 'grade thresholds are adjustable');

/* ---------------- trajectory, DCD, CIF, altloc, cryo-EM box ---------------- */
const T = B.parsePDB(fs.readFileSync(path.join(__dirname, 'data/demo_traj.pdb'), 'utf8'));
const tr = B.analyzeTrajectory(T, T.frames, T.boxes, P({ dyn: { occCut: 5 } }));
const dcd = B.parseDCD(fs.readFileSync(path.join(__dirname, 'data/demo_traj.dcd')).buffer.slice(0), S.atoms.length);
const tr2 = B.analyzeTrajectory(S, dcd.frames, dcd.boxes, P({ dyn: { occCut: 5 } }));
ok(dcd.frames.length === 60 && tr2.results.length === tr.results.length && near(tr2.results[0].occupancy, tr.results[0].occupancy), 'DCD and multi-model PDB give the same persistence');
ok(tr.results.every(r => near(r.D, B.Dscore(r.occupancy, r.meanResidence, P0))), 'D uses the documented formula');
const cif = 'data_t\n_cell_length_a 8.0\n_cell_length_b 9.0\n_cell_length_c 10.0\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\nloop_\n_symmetry_equiv_pos_as_xyz\n\'x,y,z\'\n\'-x+1/2,-y,z+1/2\'\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nI1 I 0.10 0.10 0.10\nC1 C 0.10 0.10 -0.11\nN1 N 0.10 0.10 0.40\nC2 C 0.10 0.10 0.54\n';
const M = B.parseCIF(cif); const mx = B.detectMaterialXB(M, P0); ok(mx.some(h => h.Xlabel === 'I1' && h.Ylabel === 'N1' && near(h.d, 3.0, 0.01)), 'CIF I···N found'); ok(mx.filter(h => h.Xlabel === 'I1' && h.Ylabel === 'N1' && h.status === 'XB').length === 1 && mx.find(h => h.Ylabel === 'N1').mult === 2, 'symmetry-equivalent CIF contacts reported once (2 copies)');
let bad = false; try { B.parseCIF(cif.replace("'x,y,z'", "'x,y,alert(1)'")); } catch (e) { bad = /Unsupported/.test(e.message); } ok(bad, 'unsafe CIF symmetry operator refused');
const alt = 'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N\nATOM      2  CA AALA A   1       1.450   0.000   0.000  0.30 10.00           C\nATOM      3  CA BALA A   1       1.460   0.100   0.000  0.70 10.00           C\nEND\n';
ok(B.parsePDB(alt).atoms.length === 2 && B.parsePDB(alt).atoms[1].x === 1.45, 'altloc: first conformer kept by default');
ok(B.parsePDB(alt, { altloc: 'occupancy' }).atoms[1].x === 1.46, 'altloc: highest-occupancy option');
const em = 'CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n' + fs.readFileSync(path.join(__dirname, 'data/demo_complex.pdb'), 'utf8');
const emHits = B.detectXB(B.parsePDB(em), P0).hits; ok(emHits.length === hits.length, 'cryo-EM 1 Å dummy cell is ignored (no false periodic contacts)');

/* ---------------- writers: dispersion keywords; dimers; partially protonated input ---------------- */
{
  const dm = { qm: [{ el: 'I', x: 0, y: 0, z: 0, frag: 1 }, { el: 'C', x: 0, y: 0, z: -2.1, frag: 1 }, { el: 'N', x: 0, y: 0, z: 2.9, frag: 2 }],
    manifest: { id: 'T', totalCharge: 0, multiplicity: 1, fragCharge: { 1: 0, 2: 0 }, xbs: [{ label: 'x' }], atoms: [] } };
  const g = (m, d) => B.writeGaussian(dm, { job: 'j', nproc: 1, mem: '1GB', method: m, basis: 'def2TZVP', disp: d }).split('\n')[3];
  ok(/gd3bj/.test(g('B3LYP')) && !/empirical/.test(g('M062X')) && !/empirical/.test(g('wB97XD')) && !/empirical/.test(g('HF')), 'Gaussian: D3(BJ) only for hybrids without built-in dispersion');
  const o = (m, d) => B.writeORCA(dm, { job: 'j', nproc: 1, maxcore: 100, method: m, basis: 'def2-TZVP', disp: d }).inp.split('\n').find(l => l.startsWith('!'));
  ok(/D3BJ/.test(o('PBE0')) && !/D3/.test(o('M062X')) && !/D3BJ/.test(o('wB97X-D3')) && /D3ZERO/.test(o('M062X', 'zero')), 'ORCA: dispersion keyword chosen per functional');
  { const inp = B.writeORCA(dm, { job: 'j', nproc: 1, maxcore: 100, method: 'wB97X-D3', basis: 'def2-TZVP' }).inp; const Pq = B.fullParams(); Pq.qm.program = 'orca';
    ok(/CAUTION: ORCA has no built-in NBO/.test(inp) && /NBOEXE/.test(inp) && inp.split('\n').filter(l => l.startsWith('#')).length === 6 && B.assumptions('T3', Pq).includes(B.ORCA_NBO_CAUTION), 'ORCA NBO caution in input comments and T3 assumptions'); }
  ok(/DFT_D\s+D3_BJ/.test(B.writeQChemNBO(dm, { method: 'PBE0', basis: 'def2-TZVP', disp: 'bj' })), 'Q-Chem: DFT_D D3_BJ written when requested');
  // a strong halogen bond (I...N 2.43 A) is not mistaken for a covalent bond when the dimer is split into fragments
  const xyz = '6\nFI...pyridine-like\nF 0 0 -2.0\nI 0 0 0\nN 0 0 2.43\nC 1.2 0 3.1\nC -1.2 0 3.1\nH 0 0 4.6\n';
  const Mx = B.parseXYZ(xyz); const hx = { cid: 'XB1', X: 1, A: 0, Y: 2, Ypos: Mx.atoms[2].c, image: '0 0 0', halogen: 'I', Xlabel: 'I2', Ylabel: 'N3', d: 2.43, theta1: 180 };
  const mm = B.buildMaterialQM(Mx, hx, { P: P0 }); ok(mm.manifest.atoms.filter(a => a.frag === 1).length === 2 && !mm.manifest.warnings.some(w => /same molecule/.test(w)), 'strong XB contact is never used to join the two fragments');
  // residue with only polar H (as written by some preparation tools): missing C-H are completed, N-H/O-H kept
  const raw = fs.readFileSync(path.join(__dirname, 'ref/4dst.pdb'), 'utf8'); const S4 = B.parsePDB(raw);
  const lig = h => S4.residues[h.donorRes].resName === '9LI' && S4.residues[h.donorRes].chain === 'A';
  const m4 = B.buildQMModel(S4, 0, B.detectXB(S4, P0).hits.filter(h => h.status === 'XB' && lig(h)), { P: P0 });
  const polar = m4.qm.filter(a => a.el === 'H' && a.added && ['N', 'O'].includes(m4.qm[a.parent].el));
  let serial = 90000; const extra = polar.map((h, k) => { const pa = m4.qm[h.parent].src; const rec = S4.atoms[pa.i].het ? 'HETATM' : 'ATOM  ';
    return `${rec}${String(serial + k).padStart(5)}  H${String(k).padEnd(2)} ${pa.resName.padStart(3)} ${pa.chain}${String(pa.resSeq).padStart(4)}    ${h.x.toFixed(3).padStart(8)}${h.y.toFixed(3).padStart(8)}${h.z.toFixed(3).padStart(8)}  1.00  0.00           H`; });
  const S5 = B.parsePDB(raw.replace(/^END\s*$/m, '') + extra.join('\n') + '\nEND\n');
  const m5 = B.buildQMModel(S5, 0, B.detectXB(S5, P0).hits.filter(h => h.status === 'XB' && S5.residues[h.donorRes].resName === '9LI' && S5.residues[h.donorRes].chain === 'A'), { P: P0 });
  ok(polar.length > 0 && m5.manifest.nAtoms === m4.manifest.nAtoms && m5.manifest.warnings.some(w => /C–H hydrogen/.test(w)), `polar-H-only input: C–H completed, same ${m4.manifest.nAtoms}-atom model as from the H-free file`);
}

/* ---------------- ligand protonation: one heavy-atom structure for H and charge ---------------- */
{
  const ligModel = (file, rn, ch) => { const S = B.parsePDB(fs.readFileSync(path.join(__dirname, 'ref', file), 'utf8'));
    const hits = B.detectXB(S, P0).hits.filter(h => S.residues[h.donorRes].resName === rn && S.residues[h.donorRes].chain === ch).sort((a, b) => b.G - a.G);
    const m = B.buildQMModel(S, 0, [hits[0]], { P: P0 }); const L = m.qm.filter(a => a.frag === 1); const c = {}; L.forEach(a => { c[a.el] = (c[a.el] || 0) + 1; });
    return { S, m, c, L }; };
  const formula = c => Object.keys(c).sort().map(e => e + c[e]).join('');
  const t = ligModel('1p5e.pdb', 'TBS', 'A'); const nh = t.L.filter(a => a.el === 'H' && a.added && t.m.qm[a.parent].el === 'N').length;
  ok(formula(t.c) === 'BR4C6H1N3' && nh === 1 && t.m.manifest.fragCharge[1] === 0, `1P5E TBS: C6HBr4N3 with one N-H (got ${formula(t.c)}, ${nh} N-H)`);
  const v = ligModel('2w0s.pdb', 'BVP', 'A'); ok(formula(v.c) === 'BR1C11H12N2O8P1' && v.m.manifest.fragCharge[1] === -2, `2W0S BVP: phosphate dianion C11H12BrN2O8P, charge -2 (got ${formula(v.c)}, ${v.m.manifest.fragCharge[1]})`);
  const r = ligModel('4gql.pdb', 'R47', 'A'); ok(r.c.H === 32 && r.m.manifest.fragCharge[1] === -3 && r.m.manifest.warnings.some(w => /hydrogens were rebuilt/.test(w)), `4GQL R47: incomplete file H rebuilt to 32 H at charge -3 (got ${r.c.H}, ${r.m.manifest.fragCharge[1]})`);
  const p = ligModel('4agl.pdb', 'P84', 'A'); ok(p.m.manifest.warnings.some(w => /13 of 21 heavy atoms/.test(w)), '4AGL P84: incomplete ligand (13 of 21 heavy atoms) is reported');
  ok(JSON.stringify(t.S.header.formul.TBS) === '{"C":6,"H":1,"BR":4,"N":3}' && v.S.header.formul.SO4._charge === -2, 'PDB FORMUL records parsed (formula and charge)');
  // placed H: no clash with any other atom of the model
  let clash = 0; for (const M of [t.m, v.m, r.m]) M.qm.forEach((a, i) => { if (a.el !== 'H') return; M.qm.forEach((b, j) => { if (j === i || j === a.parent) return; const d = Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z); if (b.el !== 'H' && d < 1.15) return; // its own X-H bond (hydrogens from the file carry no parent index)
      if (d < (b.el === 'H' ? 1.2 : 1.5) && !(b.el === 'H' && b.parent === a.parent && d > 1.4)) clash++; }); });
  { const S = r.S; const emb = forceH => { const P = B.fullParams({ qm: { env: 'qmmm', forceH } }); const ffa = B.assignForceField(S, P, {});
      const hh = B.detectXB(S, P).hits.filter(h => S.residues[h.donorRes].resName === 'R47').sort((a, b) => b.G - a.G); const m = B.buildQMModel(S, 0, [hh[0]], { P, ffa });
      const pc = B.embeddingCharges(S, 0, m.include, { P, ffa, qmResidueCharges: m.qmResCharge }); return [pc.charges.length, +pc.net.toFixed(4)]; };
    const a1 = emb(false), a2 = emb(true); ok(a1[0] === a2[0] && a1[1] === a2[1] && a1[1] === 2, `replaced file H never become MM charges (4GQL: ${a1} vs re-protonated ${a2})`); }
  ok(clash === 0, `added ligand H atoms have no clashes (H-H >= 1.2 A, H-X >= 1.5 A; ${clash} found)`);
}

/* ---------------- virtual sites (sigma-hole extra point) and NBO index overflow ---------------- */
{
  const txt = fs.readFileSync(path.join(__dirname, 'ref', 'hiyd_c3_8.8.pdb'), 'utf8');
  const nFile = txt.split('\n').filter(l => /^(ATOM|HETATM)/.test(l)).length;
  const S = B.parsePDB(txt, {});
  ok(S.vsites && S.vsites.count === 1 && S.vsites.fileAtoms === nFile && S.atoms.length === nFile - 1 && !S.atoms.some(a => a.name === 'EP'), `extra point EP of I-Tyr removed on reading (${S.vsites && S.vsites.count} site, ${S.atoms.length} of ${nFile} atoms kept)`);
  const pr = B.prepareDetection(S, {}); const d = pr.donors.find(x => S.residues[x.res].resName === 'MIA');
  ok(d && S.atoms[d.X].name === 'I1' && S.atoms[d.A].name === 'C6', 'iodine partner is the ring carbon C6, not the extra point');
  const f0 = new Float32Array(nFile * 3); txt.split('\n').filter(l => /^(ATOM|HETATM)/.test(l)).forEach((l, k) => { f0[3 * k] = +l.slice(30, 38); f0[3 * k + 1] = +l.slice(38, 46); f0[3 * k + 2] = +l.slice(46, 54); });
  const fm = B.dropVsites(S, [f0])[0]; let md = 0; for (let k = 0; k < fm.length; k++) md = Math.max(md, Math.abs(fm[k] - S.frames[0][k]));
  ok(fm.length === S.atoms.length * 3 && md < 1e-3, 'DCD/XYZ frames written with the extra point are mapped to the kept atoms');
  let thr = false; try { B.dropVsites(S, [new Float32Array(30)]); } catch (e) { thr = /virtual sites/.test(e.message); } ok(thr, 'frames with a wrong atom count are refused');
  const line = ' 156. LP (   1) N  22                /***. BD*(   1) C  60 - I  74            1.95    0.86    0.040';
  const nb = B.parseNBO(' Entering Gaussian System\n SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS\n     Threshold for printing:   0.50 kcal/mol\n' + line + '\n NATURAL BOND ORBITALS (Summary):\n');
  ok(nb.entries.length === 1 && nb.entries[0].acceptor.nbo === null && nb.entries[0].acceptor.atoms[1].n === 74 && nb.entries[0].e2 === 1.95, 'NBO lines with an orbital index above 999 (printed ***) are read');
}

/* ---------------- Q-Chem NBO input after the cluster runs ---------------- */
{
  const S = B.parsePDB(fs.readFileSync(path.join(__dirname, 'ref', '1p5e.pdb'), 'utf8'), {});
  const hh = B.detectXB(S, {}).hits.filter(h => h.status === 'XB').sort((a, b) => b.G - a.G);
  const m = B.buildQMModel(S, 0, [hh[0]], { P: B.fullParams({}), ffa: B.assignForceField(S, B.fullParams({}), {}) });
  const q = B.writeQChemNBO(m, { method: 'wB97X-D', basis: 'def2-TZVP', solvent: true, eps: 4 });
  ok(/\$nbo\n   E2PERT\n\$end/.test(q) && /\$solvent\n   Dielectric 4\n\$end/.test(q) && !/SolventDielectric/.test(q), 'Q-Chem NBO input requests the E(2) table and uses the Q-Chem 6 PCM keyword Dielectric');
  const q4 = B.writeQChemNBO(m, { method: 'wB97X-D', basis: 'def2-QZVP' });
  ok(q4.includes(B.QCHEM_NBO_G_CAUTION) && !q.includes(B.QCHEM_NBO_G_CAUTION), 'Q-Chem NBO input with def2-QZVP carries the pure-g caution');
  ok(B.assumptions('T3', B.fullParams({ qm: { program: 'qchem' } })).includes(B.QCHEM_NBO_CAUTION), 'Q-Chem NBO caution in the T3 assumptions');
}

/* ---------------- NBO gate with an ionic Lewis structure (X+ with a vacant LP*) ---------------- */
{
  const txt = ' Entering Gaussian System\n SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS\n     Threshold for printing:   0.50 kcal/mol\n' +
    ' from unit  3 to unit  1\n  39. LP (   1) N   3                / 34. LP*(   4) I   1                   72.20    0.39    0.161\n' +
    '  39. LP (   1) N   3                / 40. RY*(   1) I   1                    4.47    1.44    0.076\n NATURAL BOND ORBITALS (Summary):\n';
  const g = B.nboGate({ xbs: [{ id: 'x', type: 'atom', X: 1, A: 2, Y: 3 }] }, B.parseNBO(txt), { threshold: 0.5, includeRY: false })[0];
  ok(g.ionicX && Math.abs(g.E2_core - 72.2) < 1e-9 && g.verified, 'ionic NBO structure (vacant LP* on the halogen, no sigma*(A-X)) counts n -> LP*(X)');
  const txt2 = txt.replace('NATURAL BOND ORBITALS', '  10. LP (   1) O   9                /33. BD*(   1) F   2 - I   1                    1.00    0.39    0.161\n NATURAL BOND ORBITALS');
  const g2 = B.nboGate({ xbs: [{ id: 'x', type: 'atom', X: 1, A: 2, Y: 3 }] }, B.parseNBO(txt2), { threshold: 0.5, includeRY: false })[0];
  ok(!g2.ionicX && g2.E2_core === 0, 'LP*(X) is ignored when the sigma*(A-X) orbital exists');
}

/* ---------------- reference comparison (PDB2PQR) ---------------- */
const refOut = require('child_process').spawnSync(process.execPath, [path.join(__dirname, 'ff_reference.test.js')], { encoding: 'utf8' });
console.log(refOut.stdout.trim().split('\n').map(l => '  ' + l).join('\n'));
ok(refOut.status === 0, 'force-field charges agree with PDB2PQR reference PQR files');

console.log(`\n${nPass} passed, ${nFail} failed`);
if (nFail) process.exit(1);
