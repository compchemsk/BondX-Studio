/* =====================================================================
 * BondX Engine v2.0  —  Tiered halogen-bond (XB) characterisation core
 * Authors: Soumyajit Karmakar and Sabyashachi Mishra (IIT Kharagpur)
 * Pure functions, no DOM. Runs in the browser and in Node (for testing).
 *
 * Tiers:
 *   T0  Structure input (PDB / PQR / multi-model)
 *   T1  Geometric detection (static PDB / batch / CIF)             -> G
 *   T2  Dynamic persistence over MD trajectories                  -> D
 *   T3a QM model: force-field protonation + charges, QM cluster or QM/MM
 *   T3b NBO E(2) n->sigma* gate                                     -> E
 *   T4  ALMO-EDA (Q-Chem EDA2) POL vs CT decomposition             -> C
 *   T5  Fusion: BondX Quality Index (BXQI) + mechanistic class
 *
 * Every number the engine uses lives in DEFAULTS and can be overridden.
 * PARAM_SPEC documents each one (label, unit, meaning); assumptions(stage, P)
 * states in words what the current settings assume.
 * Force-field charges come from BondXFF (ff_data.js), generated from the
 * OpenMM XML distributions of AMBER ff14SB/ff19SB/OL15/OL3 and CHARMM36.
 * ===================================================================== */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory(require('./ff_data.js'));
  else root.BondX = factory(root.BondXFF);
})(typeof self !== 'undefined' ? self : this, function (FFDATA) {
  'use strict';
  const VERSION = '2.0.0';
  const QCHEM_G_BASIS = /QZV|cc-pVQZ|cc-pV5Z/i;
  const QCHEM_NBO_G_CAUTION = 'CAUTION: the built-in NBO 5 of Q-Chem stops for basis sets with pure g functions (e.g. def2-QZVP). Use def2-TZVP(D) or Gaussian for the NBO job.';
  const ORCA_NBO_CAUTION = 'ORCA has no built-in NBO. The ORCA NBO job needs NBO 7, a licensed program, set through the NBOEXE variable. Users with NBO 7 can use this route. This route is not tested in this version of BondX. Gaussian (pop=nbo) is recommended.';
  const QCHEM_NBO_CAUTION = 'In our tests, Q-Chem 6.2 with NBO 1 printed the natural bond orbitals without a second-order perturbation table. The NBO gate needs this table. BondX now adds a $nbo section with E2PERT. This addition is not yet tested. Gaussian (pop=nbo) is the tested NBO route. Q-Chem remains the program for ALMO-EDA.';
  const G16_NBO7_CAUTION = 'pop=nbo7 needs NBO 7, a licensed program linked to Gaussian. This option is not tested in this version of BondX. The default pop=nbo uses the NBO program built into Gaussian.';

  /* ------------------------------------------------------------------ */
  /*  Element data (not user parameters)                                 */
  /* ------------------------------------------------------------------ */
  const COV = { H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57, P: 1.07, S: 1.05, CL: 1.02, BR: 1.20, I: 1.39, SE: 1.20, B: 0.84, SI: 1.11,
    NA: 1.66, K: 2.03, MG: 1.41, CA: 1.76, ZN: 1.22, FE: 1.32, CU: 1.32, MN: 1.39, CO: 1.26, NI: 1.24, HG: 1.32, CD: 1.44 };
  const METALS = new Set(['NA', 'K', 'MG', 'CA', 'ZN', 'FE', 'CU', 'MN', 'CO', 'NI', 'HG', 'CD', 'LI', 'SR', 'BA', 'PT', 'AU', 'AG', 'PB', 'AL', 'GA', 'MO', 'W', 'V', 'CR', 'RB', 'CS']);
  const Z = { H: 1, HE: 2, LI: 3, BE: 4, B: 5, C: 6, N: 7, O: 8, F: 9, NE: 10, NA: 11, MG: 12, AL: 13, SI: 14, P: 15, S: 16, CL: 17, AR: 18, K: 19, CA: 20,
    SC: 21, TI: 22, V: 23, CR: 24, MN: 25, FE: 26, CO: 27, NI: 28, CU: 29, ZN: 30, GA: 31, GE: 32, AS: 33, SE: 34, BR: 35, KR: 36, RB: 37, SR: 38,
    MO: 42, AG: 47, CD: 48, SN: 50, SB: 51, TE: 52, I: 53, XE: 54, CS: 55, BA: 56, PT: 78, AU: 79, HG: 80, PB: 82 };
  const ZSYM = {}; Object.keys(Z).forEach(k => { ZSYM[Z[k]] = k; });
  const WATER = new Set(['HOH', 'WAT', 'TIP3', 'TIP', 'SOL', 'TIP4', 'H2O', 'DOD', 'T3P', 'T4P', 'SPC']);
  const STD_AA = new Set(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX', 'CYM', 'ASH', 'GLH', 'LYN', 'MSE']);
  const NUCLEIC = new Set(['DA', 'DC', 'DG', 'DT', 'DU', 'A', 'C', 'G', 'U', 'RA', 'RC', 'RG', 'RU']);
  const AA_ALIAS = { HID: 'HIS', HIE: 'HIS', HIP: 'HIS', HSD: 'HIS', HSE: 'HIS', HSP: 'HIS', CYX: 'CYS', CYM: 'CYS', ASH: 'ASP', GLH: 'GLU', LYN: 'LYS', MSE: 'MET' };
  const XB_HALOGENS = new Set(['CL', 'BR', 'I']);
  const HAL_NAME = { CL: 'Cl', BR: 'Br', I: 'I', F: 'F' };

  /* ------------------------------------------------------------------ */
  /*  DEFAULTS: every adjustable number, in one place                    */
  /* ------------------------------------------------------------------ */
  const DEFAULTS = {
    struct: { altloc: 'first', bondTol: 0.45 },
    geom: {
      // calibrated cut-offs (Å) from Karmakar & Mishra, JCIM 2026; F only used when includeF is on
      cutoffs: { I: { N: 3.97, O: 3.82, S: 4.15, AR: 4.00 }, BR: { N: 3.844, O: 3.6812, S: 4.02, AR: 3.87 }, CL: { N: 3.772, O: 3.6092, S: 3.948, AR: 3.798 }, F: { N: 3.32, O: 3.29, S: 3.57, AR: 3.40 } },
      theta1Min: 140, theta1Max: 180, theta2Use: false, theta2Min: 90, theta2Max: 150, piNormalMax: 60,
      includeF: false, includePi: true, includeHisPi: false,
      acceptorScope: { protein: true, nucleic: true, water: false, otherLigand: false },
      excludeNonLewis: false, reportNearMiss: false, selectedChains: null
    },
    score: { vdw: { H: 1.20, C: 1.70, N: 1.55, O: 1.52, F: 1.47, P: 1.80, S: 1.80, CL: 1.75, BR: 1.85, I: 1.98, SE: 1.90, B: 1.92, SI: 2.10 }, ringRadius: 1.77, Rbest: 0.85, th2Opt: 120, th2Width: 40 },
    dyn: { stride: 1, occCut: 10, dt: 1, tauScale: 5, dExp: 0.25 },
    mat: { images: 1, minDist: 1.6 },
    ff: { name: 'ff14SB', water: 'TIP3P', hisDefault: 'HIE', termini: 'charged', chainBreaks: 'internal', ssCut: 2.5,
      chargeSource: 'auto', foldH: true, missing: 'redistribute', overrides: {} },
    qm: { mode: 'full', shell: 0, shellWater: false, shellIons: false, forceH: false, linkScheme: 'fixed',
      xh: { C: 1.09, N: 1.01, O: 0.96, S: 1.34, P: 1.42, B: 1.19, SE: 1.46 }, ligandCharge: 'auto',
      env: 'cluster', pcm: true, eps: 4, program: 'gaussian', method: 'wB97XD', basis: 'def2TZVP', nproc: 16, mem: '32GB', nbo7: false },
    mm: { radius: 15, selection: 'residue', boundary: 'RC', neutralizeCut: true },
    nbo: { threshold: 0.5, geomTol: 0.05, includeRY: true },
    eda: { method: 'wB97X-V', basis: 'def2-TZVPD', preset: 1, bsse: true, covp: false, pcm: false, eps: 4, memMB: 32000, chiHigh: 0.6, chiLow: 0.4 },
    bxqi: { weights: { G: 0.25, D: 0.20, E: 0.30, C: 0.25 }, E0: 2.0, C0: 10.0, gradeA: 75, gradeB: 50, gradeC: 25, requireGateForA: true,
      consE2hi: 4, consCTlo: 0.5, consE2lo: 1, consCThi: 4 }
  };
  const clone = o => JSON.parse(JSON.stringify(o));
  function mergeParams(base, over) { // deep merge, arrays and null replaced
    const out = clone(base);
    (function m(t, s) { for (const k of Object.keys(s || {})) { const v = s[k]; if (v && typeof v === 'object' && !Array.isArray(v) && t[k] && typeof t[k] === 'object' && !Array.isArray(t[k])) m(t[k], v); else t[k] = clone(v === undefined ? null : v); } })(out, over);
    return out;
  }
  function defaults() { return clone(DEFAULTS); }
  function getPath(o, p) { return p.split('.').reduce((a, k) => (a == null ? a : a[k]), o); }
  function setPath(o, p, v) { const ks = p.split('.'); let t = o; ks.slice(0, -1).forEach(k => { if (t[k] == null || typeof t[k] !== 'object') t[k] = {}; t = t[k]; }); t[ks[ks.length - 1]] = v; }

  // Parameter documentation (drives the UI forms and the parameter table)
  const S_ = (path, stage, label, unit, type, desc, extra) => Object.assign({ path, stage, label, unit, type, desc }, extra || {});
  const PARAM_SPEC = [
    S_('struct.altloc', 'T0', 'Alternate locations', '', 'select', 'Which alternate conformer to keep when a residue has several.', { options: [['first', 'first listed (A / 1 / blank)'], ['occupancy', 'highest total occupancy']] }),
    S_('struct.bondTol', 'T0', 'Bond tolerance', 'Å', 'number', 'Two atoms are bonded when d < r_cov(1) + r_cov(2) + tolerance.', { step: 0.05 }),
    S_('geom.cutoffs', 'T1', 'X···Y distance cut-offs', 'Å', 'cutoffs', 'Maximum X···Y distance per halogen and acceptor type (π = ring centroid).'),
    S_('geom.theta1Min', 'T1', 'θ₁ minimum (A–X···Y)', '°', 'number', 'Lower limit of the σ-hole angle.'),
    S_('geom.theta1Max', 'T1', 'θ₁ maximum', '°', 'number', 'Upper limit of the σ-hole angle.'),
    S_('geom.theta2Use', 'T1', 'Filter on θ₂ (X···Y–R)', '', 'bool', 'Reject contacts whose θ₂ is outside the window below.'),
    S_('geom.theta2Min', 'T1', 'θ₂ minimum', '°', 'number', ''), S_('geom.theta2Max', 'T1', 'θ₂ maximum', '°', 'number', ''),
    S_('geom.piNormalMax', 'T1', 'π acceptor: max angle to ring normal', '°', 'number', 'Largest angle α between X→centroid and the ring normal.'),
    S_('geom.includePi', 'T1', 'Aromatic π acceptors (Phe/Tyr/Trp)', '', 'bool', ''), S_('geom.includeHisPi', 'T1', 'His ring as π acceptor', '', 'bool', ''),
    S_('geom.acceptorScope.protein', 'T1', 'Protein N/O/S acceptors', '', 'bool', ''), S_('geom.acceptorScope.nucleic', 'T1', 'Nucleic-acid acceptors', '', 'bool', ''),
    S_('geom.acceptorScope.water', 'T1', 'Water O acceptors', '', 'bool', ''), S_('geom.acceptorScope.otherLigand', 'T1', 'Other ligands / cofactors as acceptors', '', 'bool', ''),
    S_('geom.excludeNonLewis', 'T1', 'Exclude poor Lewis bases', '', 'bool', 'Drop amide N–H, Lys NZ, Arg N, Trp NE1 instead of only flagging them.'),
    S_('geom.reportNearMiss', 'T1', 'Report near-misses', '', 'bool', 'Also list contacts that pass the distance but fail the angle test.'),
    S_('geom.includeF', 'T1', 'Include C–F donors', '', 'bool', 'Fluorine rarely has a σ-hole; off by default.'),
    S_('score.vdw', 'T1', 'van der Waals radii (Bondi)', 'Å', 'vdw', 'Used for R = d / (r_X + r_Y).'),
    S_('score.ringRadius', 'T1', 'Ring pseudo-radius', 'Å', 'number', 'Radius used for an aromatic centroid in R.', { step: 0.01 }),
    S_('score.Rbest', 'T1', 'R for full distance score', '', 'number', 'g_d = 1 at R ≤ this value.', { step: 0.01 }),
    S_('score.th2Opt', 'T1', 'Optimal θ₂', '°', 'number', 'Centre of the θ₂ Gaussian (sp² lone pair ≈ 120°).'),
    S_('score.th2Width', 'T1', 'θ₂ width', '°', 'number', 'Width of the θ₂ Gaussian.'),
    S_('dyn.stride', 'T2', 'Frame stride', '', 'number', 'Analyse every n-th frame.'), S_('dyn.occCut', 'T2', 'Occupancy cut-off', '%', 'number', 'Report contacts present in at least this share of analysed frames.'),
    S_('dyn.dt', 'T2', 'Time per frame', 'ps', 'number', 'Only used to convert frames to time.', { step: 'any' }),
    S_('dyn.tauScale', 'T2', 'Residence scale τ₀', 'frames', 'number', 'D = f (1 − e^(−τ/τ₀))^p.'), S_('dyn.dExp', 'T2', 'Damping exponent p', '', 'number', '', { step: 0.05 }),
    S_('mat.images', 'T1', 'Periodic image shells (CIF)', '', 'number', 'Neighbour cells searched in each direction.'), S_('mat.minDist', 'T1', 'Minimum X···Y distance (CIF)', 'Å', 'number', 'Shorter contacts are treated as bonds or disorder.', { step: 0.1 }),
    S_('ff.name', 'T3', 'Force field (charges, protonation templates)', '', 'select', 'Source of residue templates and partial charges.', { options: [['ff14SB', 'AMBER ff14SB + OL15/OL3'], ['ff19SB', 'AMBER ff19SB + OL15/OL3'], ['charmm36', 'CHARMM36']] }),
    S_('ff.water', 'T3', 'Water model', '', 'select', 'Charges of water molecules in the MM region.', { optionsFn: P => Object.keys((FFDATA[P.ff.name] || {}).waters || {}).map(k => [k, k]) }),
    S_('ff.hisDefault', 'T3', 'Histidine default', '', 'select', 'Used when the file does not fix the tautomer (no HD1/HE2 hydrogens, residue named HIS).', { options: [['HIE', 'HIE / HSE (Nε–H, neutral)'], ['HID', 'HID / HSD (Nδ–H, neutral)'], ['HIP', 'HIP / HSP (both, +1)']] }),
    S_('ff.termini', 'T3', 'Chain termini', '', 'select', 'Charge state of the first and last residue of each protein chain.', { options: [['charged', 'charged (NH₃⁺ / COO⁻)'], ['neutral', 'uncharged (internal templates)']] }),
    S_('ff.chainBreaks', 'T3', 'Chain breaks (missing residues)', '', 'select', 'Residues next to a gap inside a chain.', { options: [['internal', 'treat as internal (no extra charge)'], ['charged', 'charged termini']] }),
    S_('ff.ssCut', 'T3', 'Disulfide S–S cut-off', 'Å', 'number', 'Cys pairs closer than this become CYX.', { step: 0.1 }),
    S_('ff.chargeSource', 'T3', 'MM charge source', '', 'select', '', { options: [['auto', 'PQR charges from a loaded PQR file, otherwise force field'], ['ff', 'force field always'], ['pqr', 'PQR file only']] }),
    S_('ff.foldH', 'T3', 'Missing hydrogens: fold charge onto parent atom', '', 'bool', 'If the file has no H, each H charge is added to its heavy atom, so residue charges stay exact.'),
    S_('ff.missing', 'T3', 'Missing heavy atoms', '', 'select', '', { options: [['redistribute', 'spread their charge over the residue (keeps net charge)'], ['ignore', 'ignore (net charge changes)']] }),
    S_('qm.mode', 'T3', 'Acceptor residue in QM region', '', 'select', '', { options: [['full', 'whole residue'], ['sidechain', 'side chain + Cα (backbone N, C, O removed; Cα capped with H)']] }),
    S_('qm.shell', 'T3', 'Extra QM shell around X', 'Å', 'number', 'Whole residues with any atom this close to a halogen join fragment 2 (0 = off).', { step: 0.5 }),
    S_('qm.shellWater', 'T3', 'Shell includes waters', '', 'bool', ''), S_('qm.shellIons', 'T3', 'Shell includes ions', '', 'bool', ''),
    S_('qm.forceH', 'T3', 'Re-protonate files with hydrogens', '', 'bool', ''),
    S_('qm.linkScheme', 'T3', 'Link-atom position', '', 'select', '', { options: [['fixed', 'fixed X–H bond length along the cut bond'], ['scaled', 'scaled: r_H = r_Q + g (r_M − r_Q), g = r₀(Q–H)/r₀(Q–M)']] }),
    S_('qm.xh', 'T3', 'X–H bond lengths', 'Å', 'xh', 'Used for added and link hydrogens.'),
    S_('qm.ligandCharge', 'T3', 'Ligand net charge', '', 'select', '', { options: [['auto', 'MOL2 total from a loaded file, otherwise the ligand structure with pH-7 rules'], ['rules', 'heavy-atom structure + pH-7 group rules (COO⁻, phosphate, sulfonate/sulfate, R₄N⁺)'], ['zero', 'zero']] }),
    S_('qm.env', 'T3', 'Environment', '', 'select', 'How the protein around the QM region is represented.', { options: [['cluster', 'QM cluster (optionally in a dielectric)'], ['qmmm', 'QM/MM single point (electrostatic embedding)']] }),
    S_('qm.pcm', 'T3', 'Cluster: implicit dielectric', '', 'bool', 'CPCM (Gaussian/ORCA) or PCM (Q-Chem) around the cluster.'),
    S_('qm.eps', 'T3', 'Dielectric constant ε', '', 'number', 'Typical protein interior: 4.', { step: 0.5 }),
    S_('qm.program', 'T3', 'QM program (NBO job)', '', 'select', '', { options: [['gaussian', 'Gaussian 16'], ['orca', 'ORCA 5/6 (needs licensed NBO 7; untested)'], ['qchem', 'Q-Chem 6']] }),
    S_('qm.method', 'T3', 'Method', '', 'text', 'Exactly as the program expects it.'), S_('qm.basis', 'T3', 'Basis set', '', 'text', 'def2 bases add def2-ECP for I automatically in Q-Chem; Gaussian/ORCA handle it internally.'),
    S_('qm.nproc', 'T3', 'Cores', '', 'number', ''), S_('qm.mem', 'T3', 'Memory', '', 'text', 'e.g. 32GB'),
    S_('qm.nbo7', 'T3', 'Gaussian: pop=nbo7 (licensed NBO 7; untested)', '', 'bool', 'Off (default): Gaussian uses its built-in NBO (pop=nbo). On: Gaussian calls NBO 7, a licensed program that must be linked to Gaussian. This option is not tested in this version of BondX.'),
    S_('mm.radius', 'T3', 'Embedding radius', 'Å', 'number', 'MM point charges within this distance of any QM atom.'),
    S_('mm.selection', 'T3', 'Embedding selection', '', 'select', '', { options: [['residue', 'whole residues (keeps residue charges integral)'], ['atom', 'individual atoms']] }),
    S_('mm.boundary', 'T3', 'QM/MM boundary charges', '', 'select', 'M1 = MM atom bonded to the QM region.', { options: [['RC', 'RC: move q(M1) to M1–M2 bond midpoints (Lin & Truhlar)'], ['Z1', 'Z1: delete q(M1)'], ['Z2', 'Z2: delete q(M1), q(M2)'], ['Z3', 'Z3: delete q(M1)–q(M3)'], ['none', 'keep all charges']] }),
    S_('mm.neutralizeCut', 'T3', 'Integer charge for cut residues', '', 'bool', 'When a residue is split between QM and MM, shift its MM atoms so the MM part carries an integer charge.'),
    S_('nbo.threshold', 'T3b', 'NBO gate threshold ΣE(2)', 'kcal/mol', 'number', 'Verified if ΣE(2) ≥ threshold.', { step: 0.05 }),
    S_('nbo.includeRY', 'T3b', 'Count n→RY*(X) terms', '', 'bool', 'Include donation into halogen Rydberg orbitals.'),
    S_('nbo.geomTol', 'T3b', 'Geometry check tolerance', 'Å', 'number', 'Largest allowed atom displacement between output and manifest.', { step: 0.01 }),
    S_('eda.method', 'T4', 'EDA method', '', 'text', 'Q-Chem METHOD keyword.'), S_('eda.basis', 'T4', 'EDA basis', '', 'text', ''),
    S_('eda.preset', 'T4', 'EDA2 preset', '', 'select', '', { options: [['1', '1 · frozen decomposition + nMDQ-uFERF POL (standard)'], ['2', '2 · frozen decomposition + AO-block ALMO POL']] }),
    S_('eda.bsse', 'T4', 'EDA_BSSE (counterpoise for CT)', '', 'bool', ''), S_('eda.covp', 'T4', 'COVP analysis', '', 'bool', ''),
    S_('eda.pcm', 'T4', 'PCM solvent', '', 'bool', ''), S_('eda.eps', 'T4', 'PCM ε', '', 'number', '', { step: 0.5 }), S_('eda.memMB', 'T4', 'MEM_TOTAL', 'MB', 'number', ''),
    S_('eda.chiHigh', 'T4', 'χ_CT for CT-driven', '', 'number', '', { step: 0.05 }), S_('eda.chiLow', 'T4', 'χ_CT for polarization-driven', '', 'number', '', { step: 0.05 }),
    S_('bxqi.weights', 'T5', 'Tier weights', '', 'weights', 'Weights of G, D, E, C in the geometric mean.'),
    S_('bxqi.E0', 'T5', 'E₀ (NBO saturation)', 'kcal/mol', 'number', 's_E = 1 − e^(−ΣE(2)/E₀)', { step: 0.5 }), S_('bxqi.C0', 'T5', 'C₀ (EDA saturation)', 'kJ/mol', 'number', 's_C = 1 − e^(−|ΔE_int|/C₀)'),
    S_('bxqi.gradeA', 'T5', 'Grade A from', '', 'number', ''), S_('bxqi.gradeB', 'T5', 'Grade B from', '', 'number', ''), S_('bxqi.gradeC', 'T5', 'Grade C from', '', 'number', ''),
    S_('bxqi.requireGateForA', 'T5', 'Grade A needs a passed NBO gate', '', 'bool', ''),
    S_('bxqi.consE2hi', 'T5', 'Consistency: large E(2)', 'kJ/mol', 'number', 'NBO E(2) is converted from kcal/mol (×4.184) before comparison.'), S_('bxqi.consCTlo', 'T5', 'Consistency: negligible CT', 'kJ/mol', 'number', '', { step: 0.1 }),
    S_('bxqi.consE2lo', 'T5', 'Consistency: small E(2)', 'kJ/mol', 'number', 'NBO E(2) converted to kJ/mol.'), S_('bxqi.consCThi', 'T5', 'Consistency: large CT', 'kJ/mol', 'number', '')
  ];

  /* ------------------------------------------------------------------ */
  /*  Assumptions in plain words (shown on every page, written to runs) */
  /* ------------------------------------------------------------------ */
  function assumptions(stage, P) {
    P = P || DEFAULTS; const g = P.geom, sc = P.score, f = P.ff, q = P.qm, m = P.mm; const ffl = (FFDATA[f.name] || {}).label || f.name; const A = [];
    if (stage === 'T0') {
      A.push('The first MODEL of a multi-model file defines the atoms. Later models serve as trajectory frames.');
      A.push(P.struct.altloc === 'occupancy' ? 'For alternate locations, the conformer with the highest total occupancy is kept.' : 'For alternate locations, the first conformer (A, 1 or blank) is kept.');
      A.push('Virtual sites are removed on reading. Examples are the σ-hole extra point (EP) of a halogen, LP sites and TIP4P M sites. DCD and XYZ frames with these sites are mapped onto the remaining atoms.');
      A.push(`Bonds come from distances: d < r_cov(1) + r_cov(2) + ${P.struct.bondTol} Å. Links between residues are limited to peptide, disulfide and phosphodiester bonds.`);
      A.push('Element symbols come from columns 77–78. Atom names supply missing element symbols.');
      A.push('Residues are classed as protein, nucleic acid, water, ion (single atom) or ligand (all other residues).');
    } else if (stage === 'T1') {
      A.push(`Donors: C-bound ${g.includeF ? 'F, ' : ''}Cl, Br and I atoms of ligands and modified residues.`);
      A.push(`A halogen bond needs d(X···Y) ≤ the pair cut-off and ${g.theta1Min}° ≤ θ₁ ≤ ${g.theta1Max}°${g.theta2Use ? ` and ${g.theta2Min}° ≤ θ₂ ≤ ${g.theta2Max}°` : ''}.`);
      A.push(`Acceptors: ${[g.acceptorScope.protein && 'protein N/O/S', g.acceptorScope.nucleic && 'nucleic-acid N/O', g.acceptorScope.water && 'water O', g.acceptorScope.otherLigand && 'other ligands', g.includePi && `aromatic rings (Phe, Tyr, Trp${g.includeHisPi ? ', His' : ''}) with α ≤ ${g.piNormalMax}°`].filter(Boolean).join(', ')}.`);
      A.push(`R = d / (r_X + r_Y) with Bondi radii. An aromatic centroid has r = ${sc.ringRadius} Å.`);
      A.push(`G = [g_d · g_θ · (½ + ½ g₂)]^⅓. g_d = 1 at R ≤ ${sc.Rbest}. g₂ is a Gaussian around θ₂ = ${sc.th2Opt}° (width ${sc.th2Width}°).`);
      A.push(g.excludeNonLewis ? 'Poor Lewis bases (amide N–H, Lys NZ, Arg N, Trp NE1) are excluded.' : 'Poor Lewis bases are kept with the flag “weak”. The NBO gate decides their fate.');
      A.push('Hydrogens are not required. Periodic images apply only to MD boxes (orthorhombic, edges ≥ 10 Å).');
    } else if (stage === 'T2') {
      A.push('Frames must have the atoms of the topology in the same order.');
      A.push(`The T1 test runs on every ${P.dyn.stride === 1 ? '' : P.dyn.stride + '-th '}frame. Contacts in fewer than ${P.dyn.occCut}% of the frames are not reported.`);
      A.push(`Residence time τ = mean length of unbroken runs. D = f·(1 − e^(−τ/${P.dyn.tauScale}))^${P.dyn.dExp}.`);
      A.push('The representative frame has the shortest d. Ties go to the largest θ₁.');
      A.push('Orthorhombic boxes (DCD unit cell or CRYST1) use the minimum-image convention. The ligand should be imaged and centred.');
    } else if (stage === 'T3') {
      A.push(`Protonation states and partial charges come from ${ffl}. ASP and GLU are anionic. LYS and ARG are cationic. CYS and TYR are neutral. HIS is ${f.hisDefault} by default. A tautomer name in the file (HID/HIE/HIP, HSD/HSE/HSP) or its hydrogens override this default. Changes for single residues override everything.`);
      A.push(`Cys pairs with S–S < ${f.ssCut} Å are disulfides (CYX). Chain termini are ${f.termini === 'charged' ? 'charged (NH₃⁺, COO⁻)' : 'uncharged'}. Residues next to a chain break are ${f.chainBreaks === 'charged' ? 'charged termini' : 'treated as internal'}.`);
      A.push('Fragment 1 of the QM region is the whole halogenated ligand. Fragment 2 holds ' + (q.mode === 'sidechain' ? 'the acceptor side chains with Cα. Backbone N, C and O atoms leave fragment 2 except as acceptors' : 'the whole acceptor residues') + `.${q.shell > 0 ? ` Residues within ${q.shell} Å of the halogen are added.` : ''}`);
      A.push(`Missing hydrogens follow the force-field template of each residue (number of H per atom) with VSEPR geometry. Cut bonds get H caps (${q.linkScheme === 'scaled' ? 'scaled link-atom position' : 'fixed X–H length'}).`);
      A.push(`Ligands: one heavy-atom structure fixes the hydrogens and the charge. Every sp² carbon gets one double bond. The keto form is preferred. An sp² ring N without a double bond carries one H. P keeps one P=O bond. S keeps up to two S=O bonds. ${q.ligandCharge === 'zero' ? 'All acids stay neutral.' : 'At pH 7, COOH, phosphate, phosphonate, sulfonate and sulfate groups lose their O–H hydrogen. Amines stay neutral.'} Hydrogens of the file stay only with a matching number. The PDB formula (FORMUL) checks the count.`);
      A.push('Fragment charges: residue formal charges come from the force-field template. The ligand charge comes from ' + (q.ligandCharge === 'zero' ? 'zero' : q.ligandCharge === 'rules' ? 'the ligand structure above' : 'a loaded MOL2 file or from the ligand structure above') + '. Every value can be overridden.');
      if (q.env === 'qmmm') {
        A.push(`QM/MM: MM point charges within ${m.radius} Å (${m.selection === 'residue' ? 'whole residues' : 'single atoms'}) come from ${f.chargeSource === 'pqr' ? 'the PQR file' : f.chargeSource === 'ff' ? ffl : 'a loaded PQR file or from ' + ffl}. Water model: ${f.water}.`);
        A.push(f.foldH ? 'For atoms without hydrogens in the file, each hydrogen charge moves onto its heavy atom. Residue charges stay exact. Local dipoles are approximate.' : 'Hydrogen charges apply only to hydrogens in the file.');
        A.push({ RC: 'Boundary: the charge of each M1 atom moves to the midpoints of its M1–M2 bonds (redistributed charge, RC).', Z1: 'Boundary: the charges of the M1 atoms are deleted (Z1).', Z2: 'Boundary: the charges of the M1 and M2 atoms are deleted (Z2).', Z3: 'Boundary: the charges of the M1–M3 atoms are deleted (Z3).', none: 'Boundary: all MM charges are kept (no correction).' }[m.boundary]);
        A.push('The embedding is electrostatic only. QM–MM van der Waals terms are not included. They do not enter the NBO analysis of a single point.');
        A.push('Ligands and cofactors in the MM region carry MOL2 charges from a loaded file. Without a MOL2 file, a charged group (carboxylate, phosphate, sulfate, R₄N⁺) carries its formal charge on its heavy atoms. Neutral groups get zero charges. Both cases are flagged.');
      } else A.push(q.pcm ? `The QM cluster sits in an implicit dielectric ε = ${q.eps} (CPCM/PCM). The dielectric represents the protein.` : 'The QM cluster is in the gas phase.');
      if (q.program === 'orca') A.push(ORCA_NBO_CAUTION);
      if (q.program === 'qchem') A.push(QCHEM_NBO_CAUTION);
      if (q.program === 'gaussian' && q.nbo7) A.push(G16_NBO7_CAUTION);
    } else if (stage === 'T3b') {
      A.push(`A contact is verified at ΣE(2) ≥ ${P.nbo.threshold} kcal/mol. The sum covers n(Y)→σ*(A–X)${P.nbo.includeRY ? ' and n(Y)→RY*(X)' : ''} for lone-pair acceptors. It covers π(ring)→σ*(A–X) for aromatic acceptors.`);
      A.push(`The output geometry must match the manifest atom by atom (same elements, deviation < ${P.nbo.geomTol} Å). A mismatch blocks the gate.`);
      A.push('NBO prints only couplings above its own print threshold (0.50 kcal/mol by default). Open-shell α and β tables are summed.');
    } else if (stage === 'T4') {
      A.push('Fragment 1 is the halogenated ligand. Fragment 2 is the acceptor environment. Both keep the charges and multiplicities of the QM model.');
      A.push(`χ_CT = E_CT / (E_POL + E_CT). Values ≥ ${P.eda.chiHigh} mean charge-transfer driven. Values ≤ ${P.eda.chiLow} mean polarization driven. Values in between mean mixed.`);
      A.push('Energies are given in kJ/mol (1 kcal = 4.184 kJ).');
    } else if (stage === 'T5') {
      const w = P.bxqi.weights;
      A.push(`BXQI = 100 · Π s_k^(w_k/Σw) over all completed tiers. Weights: G ${w.G}, D ${w.D}, E ${w.E}, C ${w.C}.`);
      A.push(`s_E = 1 − e^(−ΣE(2)/${P.bxqi.E0}). s_C = 1 − e^(−|ΔE_int|/${P.bxqi.C0}). Grades: A ≥ ${P.bxqi.gradeA}, B ≥ ${P.bxqi.gradeB}, C ≥ ${P.bxqi.gradeC}, else D.`);
      if (P.bxqi.requireGateForA) A.push('Grade A requires a passed NBO gate. A contact with a failed gate receives grade C or D.');
      A.push(`An NBO/ALMO inconsistency is flagged in two cases: E(2) > ${P.bxqi.consE2hi} kJ/mol with −E_CT < ${P.bxqi.consCTlo} kJ/mol, and E(2) < ${P.bxqi.consE2lo} kJ/mol with −E_CT > ${P.bxqi.consCThi} kJ/mol.`);
    }
    return A;
  }

  /* ------------------------------------------------------------------ */
  /*  Vector helpers                                                     */
  /* ------------------------------------------------------------------ */
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const scl = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm = a => Math.sqrt(dot(a, a));
  const unit = a => { const n = norm(a); return n > 1e-12 ? scl(a, 1 / n) : [0, 0, 0]; };
  const dist = (a, b) => norm(sub(a, b));
  const deg = r => r * 180 / Math.PI;
  function angle(a, b, c) { const ba = sub(a, b), bc = sub(c, b); const n1 = norm(ba), n2 = norm(bc); if (n1 === 0 || n2 === 0) return 0; return deg(Math.acos(Math.max(-1, Math.min(1, dot(ba, bc) / (n1 * n2))))); }
  function anyPerp(v) { const t = Math.abs(v[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]; return unit(cross(v, t)); }
  function rotAround(v, axis, ang) { const k = unit(axis), c = Math.cos(ang), s = Math.sin(ang); return add(add(scl(v, c), scl(cross(k, v), s)), scl(k, dot(k, v) * (1 - c))); }
  function centroid(pts) { const c = [0, 0, 0]; pts.forEach(p => { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; }); return scl(c, 1 / pts.length); }
  function planeNormal(pts) { const c = centroid(pts); let n = [0, 0, 0]; for (let i = 0; i < pts.length; i++) { const p = sub(pts[i], c), q = sub(pts[(i + 1) % pts.length], c); n = add(n, cross(p, q)); } return unit(n); }

  /* ------------------------------------------------------------------ */
  /*  Element inference                                                  */
  /* ------------------------------------------------------------------ */
  function inferElement(name, resName, het) {
    const n = name.trim().toUpperCase().replace(/[^A-Z]/g, ''); if (!n) return 'X';
    if (!het && STD_AA.has(resName)) return n[0];
    const two = n.slice(0, 2);
    if (two === 'CL' || two === 'BR') return two;
    if (n === resName && ['NA', 'MG', 'ZN', 'FE', 'CA', 'MN', 'CU', 'CO', 'NI', 'HG', 'CD', 'K', 'LI', 'SE'].includes(n)) return n;
    if (two === 'SE' && resName === 'MSE') return 'SE';
    return n[0];
  }

  /* ------------------------------------------------------------------ */
  /*  T0: PDB / PQR parsing                                              */
  /* ------------------------------------------------------------------ */
  function parsePDB(text, opts = {}) {
    const altPolicy = opts.altloc || DEFAULTS.struct.altloc; const bondTol = opts.bondTol ?? DEFAULTS.struct.bondTol;
    const lines = text.split(/\r?\n/);
    // pass 1 (altloc = occupancy): choose per residue the alternate location id with the highest summed occupancy
    let keepAlt = null;
    if (altPolicy === 'occupancy') {
      const occ = new Map(); let inModel2 = false, models = 0;
      for (const L of lines) { const rec = L.slice(0, 6).trim().toUpperCase(); if (rec === 'MODEL') { models++; if (models > 1) inModel2 = true; } if (inModel2) break;
        if ((rec === 'ATOM' || rec === 'HETATM') && L.length >= 54) { const alt = L.slice(16, 17); if (alt === ' ') continue; const rk = L.slice(17, 27); const o = parseFloat(L.slice(54, 60)) || 0;
          const m = occ.get(rk) || {}; m[alt] = (m[alt] || 0) + o; occ.set(rk, m); } }
      keepAlt = new Map(); for (const [rk, m] of occ) keepAlt.set(rk, Object.entries(m).sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1))[0][0]);
    }
    const altOK = (alt, L) => { if (!alt || alt === ' ') return true; if (keepAlt) return keepAlt.get(L.slice(17, 27)) === alt; return alt === 'A' || alt === '1'; };
    const atoms = []; const frames = []; let cur = null; let modelCount = 0; let firstModelDone = false;
    const header = { title: '', compound: {}, idCode: '', resolution: null, method: '', formul: {} }; const formulTxt = {};
    const boxes = []; let box = null; let isPQR = !!opts.pqr; let molId = null; const molNames = {};
    for (const L of lines) {
      const rec = L.slice(0, 6).trim().toUpperCase();
      if (rec === 'HEADER') header.idCode = L.slice(62, 66).trim();
      else if (rec === 'TITLE') header.title += L.slice(10).trim() + ' ';
      else if (rec === 'EXPDTA') header.method = L.slice(10).trim();
      else if (rec === 'REMARK' && /RESOLUTION\.\s+([\d.]+)/.test(L)) header.resolution = parseFloat(RegExp.$1);
      else if (rec === 'COMPND') {
        const t = L.slice(10); let m;
        if ((m = /MOL_ID:\s*(\d+)/.exec(t))) molId = m[1];
        if ((m = /MOLECULE:\s*([^;]+)/.exec(t))) molNames[molId] = m[1].trim();
        if ((m = /CHAIN:\s*([^;]+)/.exec(t))) m[1].split(',').forEach(c => { header.compound[c.trim()] = molNames[molId] || ''; });
      } else if (rec === 'FORMUL') {
        const rn = L.slice(12, 15).trim(); if (rn) formulTxt[rn] = (formulTxt[rn] || '') + ' ' + L.slice(19, 70);
      } else if (rec === 'CRYST1') {
        box = [parseFloat(L.slice(6, 15)), parseFloat(L.slice(15, 24)), parseFloat(L.slice(24, 33)), parseFloat(L.slice(33, 40)), parseFloat(L.slice(40, 47)), parseFloat(L.slice(47, 54))];
      } else if (rec === 'MODEL') { modelCount++; cur = []; }
      else if (rec === 'ENDMDL' || (rec === 'END' && cur && cur.length)) {
        if (cur && cur.length) { frames.push(Float32Array.from(cur)); boxes.push(box); firstModelDone = true; }
        cur = null;
      } else if (rec === 'ATOM' || rec === 'HETATM') {
        let x, y, z, name, resName, chain, resSeq, iCode, alt, el, occ = 1, b = 0, q = null, serial;
        const toks = L.trim().split(/\s+/);
        if (isPQR || (L.length < 54 && toks.length >= 10)) {
          isPQR = true; const n = toks.length; q = parseFloat(toks[n - 2]); const r = parseFloat(toks[n - 1]);
          z = parseFloat(toks[n - 3]); y = parseFloat(toks[n - 4]); x = parseFloat(toks[n - 5]);
          serial = parseInt(toks[1]); name = toks[2]; resName = toks[3];
          if (n === 11) { chain = toks[4]; resSeq = parseInt(toks[5]); } else { chain = ' '; resSeq = parseInt(toks[4]); }
          iCode = ''; alt = ''; el = null; b = r;
        } else {
          serial = parseInt(L.slice(6, 11)); name = L.slice(12, 16); alt = L.slice(16, 17); resName = L.slice(17, 21).trim();
          chain = L.slice(21, 22); resSeq = parseInt(L.slice(22, 26)); iCode = L.slice(26, 27).trim();
          x = parseFloat(L.slice(30, 38)); y = parseFloat(L.slice(38, 46)); z = parseFloat(L.slice(46, 54));
          occ = parseFloat(L.slice(54, 60)); if (isNaN(occ)) occ = 1; b = parseFloat(L.slice(60, 66)) || 0;
          el = L.slice(76, 78).trim().toUpperCase() || null;
          if (L.length >= 80) { const cs = L.slice(78, 80).trim(); if (cs) { const sg = cs.endsWith('-') ? -1 : 1; q = sg * (parseInt(cs) || 0); } }
        }
        if (!altOK(alt, L)) continue;
        if (cur === null && modelCount === 0) cur = [];
        if (!firstModelDone) {
          const het = rec === 'HETATM';
          const element = (el && /^[A-Z]{1,2}$/.test(el)) ? el : inferElement(name, resName.toUpperCase(), het);
          atoms.push({ i: atoms.length, serial, name: name.trim(), resName: resName.toUpperCase(), chain: chain === ' ' ? '_' : chain, resSeq, iCode, el: element, het, occ, b, q, x, y, z });
        }
        if (cur) cur.push(x, y, z);
      }
    }
    if (cur && cur.length) { frames.push(Float32Array.from(cur)); boxes.push(box); }
    if (!frames.length && atoms.length) frames.push(Float32Array.from(atoms.flatMap(a => [a.x, a.y, a.z])));
    const nA = atoms.length; const good = []; const goodBoxes = [];
    frames.forEach((f, k) => { if (f.length === nA * 3) { good.push(f); goodBoxes.push(boxes[k] || null); } });
    header.title = header.title.trim();
    // FORMUL: formula of one molecule of each HET group (water and ions included), e.g. "2(C6 H BR4 N3)" -> {C:6,H:1,BR:4,N:3}
    for (const [rn, t] of Object.entries(formulTxt)) {
      const body = (/\(([^)]*)\)/.exec(t) || [null, t])[1]; const f = {};
      body.replace(/\b([A-Z]{1,2})(\d*)\b/g, (_, e, n) => { if (/^\d/.test(e)) return; f[e] = (f[e] || 0) + (n ? parseInt(n) : 1); return ''; });
      const qm = /(\d+)([+-])\s*$/.exec(body.trim()); if (qm) f._charge = (qm[2] === '-' ? -1 : 1) * parseInt(qm[1]);
      if (Object.keys(f).length) header.formul[rn] = f;
    }
    // Virtual sites (sigma-hole extra points EP, lone-pair sites LP, TIP4P/TIP4P-Ew M sites) carry charge but no electrons or nucleus.
    // They are removed here; S.vsites keeps the file atom count and the kept indices so trajectories with the sites can be mapped.
    const isV = a => /^(EP|LP)\d*$/i.test(a.name) || ['EP', 'LP', 'XX'].includes(a.el) || ((a.name === 'EPW' || a.name === 'MW') && /^(WAT|HOH|TIP4|TP4|T4E|SOL)$/.test(a.resName));
    let vsites = null;
    if (atoms.some(isV)) {
      const keepIdx = atoms.filter(a => !isV(a)).map(a => a.i); const names = [...new Set(atoms.filter(isV).map(a => a.resName + ':' + a.name))];
      vsites = { count: atoms.length - keepIdx.length, fileAtoms: atoms.length, keep: keepIdx, names };
      const sub = f => { const o = new Float32Array(keepIdx.length * 3); keepIdx.forEach((k, n) => { o[3 * n] = f[3 * k]; o[3 * n + 1] = f[3 * k + 1]; o[3 * n + 2] = f[3 * k + 2]; }); return o; };
      for (let k = 0; k < good.length; k++) good[k] = sub(good[k]);
      const kept = keepIdx.map((k, n) => Object.assign(atoms[k], { i: n })); atoms.length = 0; kept.forEach(a => atoms.push(a));
    }
    const S = { atoms, frames: good, boxes: goodBoxes, header, isPQR, hasH: atoms.some(a => a.el === 'H'), bondTol, altloc: altPolicy, vsites };
    buildTopology(S);
    return S;
  }

  /* ------------------------------------------------------------------ */
  /*  Topology: residues, classes, bonds, chain order, termini           */
  /* ------------------------------------------------------------------ */
  function resKey(a) { return a.chain + '|' + a.resSeq + '|' + (a.iCode || '') + '|' + a.resName; }
  function buildTopology(S) {
    const res = []; const map = new Map();
    for (const a of S.atoms) {
      const k = resKey(a); let r = map.get(k);
      if (!r) { r = { idx: res.length, key: k, resName: a.resName, chain: a.chain, resSeq: a.resSeq, iCode: a.iCode, atoms: [], het: a.het }; map.set(k, r); res.push(r); }
      r.atoms.push(a.i); a.res = r.idx;
    }
    for (const r of res) {
      const rn = r.resName;
      if (WATER.has(rn)) r.cls = 'water';
      else if (STD_AA.has(rn)) r.cls = 'protein';
      else if (NUCLEIC.has(rn)) r.cls = 'nucleic';
      else if (r.atoms.length === 1) r.cls = 'ion';
      else r.cls = 'ligand';
      r.halogens = r.atoms.filter(i => ['F', 'CL', 'BR', 'I'].includes(S.atoms[i].el)).map(i => S.atoms[i].el);
      r.byName = {}; r.atoms.forEach(i => { r.byName[S.atoms[i].name] = i; });
    }
    S.residues = res;
    S.bonds = perceiveBonds(S, S.frames[0], S.bondTol ?? DEFAULTS.struct.bondTol);
    // chain order and termini (protein and nucleic acids)
    const bonded = (i, j) => i !== undefined && j !== undefined && S.bonds[i].includes(j);
    const byChain = {}; res.forEach(r => { if (r.cls === 'protein' || r.cls === 'nucleic') (byChain[r.chain] = byChain[r.chain] || []).push(r); });
    for (const list of Object.values(byChain)) {
      list.forEach((r, k) => {
        const prev = list[k - 1], next = list[k + 1];
        if (r.cls === 'protein') {
          const N = r.byName.N, C = r.byName.C;
          const linkPrev = N !== undefined && S.bonds[N].some(j => S.atoms[j].res !== r.idx && S.atoms[j].name === 'C');   // peptide bond or ACE cap
          const linkNext = C !== undefined && S.bonds[C].some(j => S.atoms[j].res !== r.idx && S.atoms[j].name === 'N');   // peptide bond or NME/NHE cap
          r.nterm = linkPrev ? null : (k === 0 || !prev || prev.cls !== 'protein' ? 'start' : 'break');
          r.cterm = linkNext ? null : (r.byName.OXT !== undefined || r.byName.OT2 !== undefined || k === list.length - 1 || !next || next.cls !== 'protein' ? 'end' : 'break');
          void bonded;
        } else {
          const linkPrev = prev && prev.cls === 'nucleic' && bonded(r.byName.P, prev.byName["O3'"]);
          const linkNext = next && next.cls === 'nucleic' && bonded(r.byName["O3'"], next.byName.P);
          r.nterm = linkPrev ? null : (k === 0 ? 'start' : 'break'); r.cterm = linkNext ? null : (k === list.length - 1 ? 'end' : 'break');
        }
      });
    }
  }
  function coordOf(frame, i) { return [frame[3 * i], frame[3 * i + 1], frame[3 * i + 2]]; }
  function makeGrid(frame, idxs, cell) { const g = new Map(); for (const i of idxs) { const k = Math.floor(frame[3 * i] / cell) + ',' + Math.floor(frame[3 * i + 1] / cell) + ',' + Math.floor(frame[3 * i + 2] / cell); let b = g.get(k); if (!b) { b = []; g.set(k, b); } b.push(i); } return { g, cell }; }
  function gridQuery(grid, p, r, cb) { const c = grid.cell, n = Math.ceil(r / c); const ix = Math.floor(p[0] / c), iy = Math.floor(p[1] / c), iz = Math.floor(p[2] / c);
    for (let dx = -n; dx <= n; dx++) for (let dy = -n; dy <= n; dy++) for (let dz = -n; dz <= n; dz++) { const b = grid.g.get((ix + dx) + ',' + (iy + dy) + ',' + (iz + dz)); if (b) for (const j of b) cb(j); } }
  function perceiveBonds(S, frame, tol) {
    const nb = S.atoms.map(() => []); const grid = makeGrid(frame, S.atoms.map(a => a.i), 2.6);
    for (const a of S.atoms) {
      if (METALS.has(a.el)) continue;
      const ra = COV[a.el] || 0.77; const pa = coordOf(frame, a.i); const ca = S.residues[a.res].cls;
      gridQuery(grid, pa, 2.6, j => {
        if (j <= a.i) return; const b = S.atoms[j]; if (METALS.has(b.el)) return;
        const cb = S.residues[b.res].cls;
        if (a.res !== b.res && (ca === 'water' || cb === 'water')) return;
        if (a.el === 'H' && b.el === 'H') return;
        const d = dist(pa, coordOf(frame, j)); const rb = COV[b.el] || 0.77;
        if (d > 0.4 && d < ra + rb + tol) {
          if (a.res !== b.res && ca !== 'ligand' && cb !== 'ligand') {
            const ok = (a.name === 'C' && b.name === 'N') || (a.name === 'N' && b.name === 'C') || (a.name === 'SG' && b.name === 'SG') || (a.name === "O3'" && b.name === 'P') || (a.name === 'P' && b.name === "O3'");
            if (!ok) return;
          }
          nb[a.i].push(j); nb[j].push(a.i);
        }
      });
    }
    return nb;
  }

  function structureReport(S) {
    const chains = {};
    for (const r of S.residues) {
      const c = chains[r.chain] || (chains[r.chain] = { chain: r.chain, aa: [], waters: 0, ligands: new Set(), ions: new Set(), nucleic: 0, halogenLigands: new Set(), breaks: 0 });
      if (r.cls === 'protein') { c.aa.push(r.resSeq); if (r.nterm === 'break') c.breaks++; }
      else if (r.cls === 'water') c.waters++; else if (r.cls === 'ion') c.ions.add(r.resName); else if (r.cls === 'nucleic') c.nucleic++;
      else if (r.cls === 'ligand') { c.ligands.add(r.resName); if (r.halogens.some(h => h !== 'F')) c.halogenLigands.add(r.resName); }
    }
    return Object.values(chains).map(c => ({ chain: c.chain, protein: c.aa.length > 0, molecule: S.header.compound[c.chain] || (c.aa.length ? 'Protein' : ''),
      length: c.aa.length, range: c.aa.length ? `${c.aa[0]} – ${c.aa[c.aa.length - 1]}` : '', waters: c.waters, nucleic: c.nucleic, breaks: c.breaks,
      ligands: [...c.ligands], ions: [...c.ions], halogenLigands: [...c.halogenLigands] }));
  }
  function filterChains(S, chains) { if (!chains || !chains.length) return S; return rebuildSubset(S, S.atoms.filter(a => chains.includes(a.chain)).map(a => a.i)); }
  function rebuildSubset(S, indices) {
    const atoms = indices.map((o, n) => Object.assign({}, S.atoms[o], { i: n }));
    const frames = S.frames.map(f => { const out = new Float32Array(indices.length * 3); indices.forEach((o, n) => { out[3 * n] = f[3 * o]; out[3 * n + 1] = f[3 * o + 1]; out[3 * n + 2] = f[3 * o + 2]; }); return out; });
    const T = { atoms, frames, boxes: S.boxes, header: S.header, isPQR: S.isPQR, hasH: atoms.some(a => a.el === 'H'), bondTol: S.bondTol, altloc: S.altloc };
    buildTopology(T); return T;
  }
  function writePDB(S, frame, indices, remark) {
    frame = frame || S.frames[0]; indices = indices || S.atoms.map(a => a.i);
    const out = []; if (remark) remark.split('\n').forEach(l => out.push(('REMARK 999 ' + l).slice(0, 80)));
    let n = 1;
    for (const i of indices) {
      const a = S.atoms[i]; const p = coordOf(frame, i);
      const nm = a.name.length < 4 && a.el.length === 1 ? ' ' + a.name.padEnd(3) : a.name.padEnd(4);
      out.push((a.het ? 'HETATM' : 'ATOM  ') + String(n++ % 100000).padStart(5) + ' ' + nm + ' ' + (a.resName.length <= 3 ? a.resName.padStart(3) + ' ' : a.resName.slice(0, 4)) + (a.chain === '_' ? ' ' : a.chain) +
        String(a.resSeq).padStart(4) + (a.iCode || ' ') + '   ' + p[0].toFixed(3).padStart(8) + p[1].toFixed(3).padStart(8) + p[2].toFixed(3).padStart(8) +
        (a.occ ?? 1).toFixed(2).padStart(6) + (a.b || 0).toFixed(2).padStart(6) + '          ' + a.el.padStart(2));
    }
    out.push('END'); return out.join('\n') + '\n';
  }

  /* ------------------------------------------------------------------ */
  /*  T1: geometric halogen-bond detection                               */
  /* ------------------------------------------------------------------ */
  const LEWIS = { _BB: { O: 'strong', OXT: 'strong', N: 'weak' }, ARG: { NE: 'weak', NH1: 'weak', NH2: 'weak' }, LYS: { NZ: 'weak' }, TRP: { NE1: 'weak' },
    ASN: { ND2: 'weak', OD1: 'strong' }, GLN: { NE2: 'weak', OE1: 'strong' }, HIS: { ND1: 'strong', NE2: 'strong' }, PRO: { N: 'weak' } };
  const RINGS = { PHE: [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']], TYR: [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
    TRP: [['CD2', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3'], ['CG', 'CD1', 'NE1', 'CE2', 'CD2']], HIS: [['CG', 'ND1', 'CE1', 'NE2', 'CD2']] };
  function resLabel(S, r) { const R = S.residues[r]; return `${R.resName}${R.resSeq}${R.iCode || ''}:${R.chain}`; }
  function atomLabel(S, i) { const a = S.atoms[i]; return `${resLabel(S, a.res)}@${a.name}`; }
  function halogenPartner(S, xi) {
    const px0 = coordOf(S.frames[0], xi); const rx = COV[S.atoms[xi].el] || 1.2;
    const nb = S.bonds[xi].filter(j => S.atoms[j].el !== 'H' && S.atoms[j].res === S.atoms[xi].res)
      .sort((p, q) => dist(px0, coordOf(S.frames[0], p)) / (rx + (COV[S.atoms[p].el] || 0.8)) - dist(px0, coordOf(S.frames[0], q)) / (rx + (COV[S.atoms[q].el] || 0.8)));
    if (nb.length) return nb[0];
    const R = S.residues[S.atoms[xi].res]; let best = -1, bd = 1e9;
    for (const j of R.atoms) { if (j === xi) continue; const d = dist(px0, coordOf(S.frames[0], j)); if (d < bd) { bd = d; best = j; } }
    return best;
  }
  function lewisClass(S, yi) {
    const a = S.atoms[yi]; const rn = AA_ALIAS[a.resName] || a.resName; const R = S.residues[a.res];
    if (R.cls === 'protein') { const t = (LEWIS[rn] && LEWIS[rn][a.name]) || LEWIS._BB[a.name]; return t || 'strong'; }
    if (a.el === 'N') { const hv = S.bonds[yi].filter(j => S.atoms[j].el !== 'H').length; return hv >= 3 ? 'weak' : 'strong'; }
    return 'strong';
  }
  function acceptorClass(S, yi) {
    const a = S.atoms[yi]; const R = S.residues[a.res]; const rn = AA_ALIAS[a.resName] || a.resName;
    if (R.cls === 'water') return 'water O'; if (R.cls === 'nucleic') return 'nucleic ' + a.el;
    if (R.cls === 'protein') {
      if (a.name === 'O' || a.name === 'OXT') return 'backbone C=O'; if (a.name === 'N') return 'backbone N–H';
      if (['OD1', 'OD2', 'OE1', 'OE2'].includes(a.name)) return (rn === 'ASP' || rn === 'GLU') ? 'carboxylate O' : 'amide C=O';
      if (['OG', 'OG1', 'OH'].includes(a.name)) return 'hydroxyl O'; if (a.name === 'SD') return 'thioether S'; if (a.name === 'SG') return 'thiol S';
      if (rn === 'HIS') return 'imidazole N'; return 'side-chain ' + a.el;
    }
    return 'ligand ' + a.el;
  }
  function ringsOf(S, residueIdx, includeHis) {
    const R = S.residues[residueIdx]; const rn = AA_ALIAS[R.resName] || R.resName; const out = [];
    const defs = RINGS[rn]; if (!defs || (rn === 'HIS' && !includeHis)) return out;
    defs.forEach((names, k) => { const ids = names.map(n => R.byName[n]).filter(i => i !== undefined); if (ids.length >= names.length - 1) out.push({ ids, label: rn === 'TRP' ? (k === 0 ? 'π(6-ring)' : 'π(5-ring)') : 'π(ring)' }); });
    return out;
  }
  function geometricScore(el, accEl, d, th1, th2, P, isPi, piAng) {
    const g = P.geom, sc = P.score;
    const sum = (sc.vdw[el] || 1.9) + (isPi ? sc.ringRadius : (sc.vdw[accEl] || 1.6));
    const cut = (g.cutoffs[el] || {})[isPi ? 'AR' : accEl] || sum;
    const R = d / sum, Rcut = cut / sum, Rbest = sc.Rbest;
    const gd = Math.max(0, Math.min(1, (Rcut - R) / Math.max(1e-6, Rcut - Rbest)));
    const gt = Math.max(0, Math.min(1, (th1 - g.theta1Min) / Math.max(1e-6, 180 - g.theta1Min)));
    let g2 = 1;
    if (isPi && piAng != null) g2 = Math.max(0, Math.min(1, 1 - piAng / g.piNormalMax));
    else if (th2 != null) g2 = Math.exp(-Math.pow((th2 - sc.th2Opt) / sc.th2Width, 2));
    return { R, Rcut, sumVdw: sum, gd, gt, g2, G: (gd > 0 && gt > 0) ? Math.pow(gd * gt * (0.5 + 0.5 * g2), 1 / 3) : 0 };
  }
  function fullParams(P) { return mergeParams(DEFAULTS, P || {}); }
  function prepareDetection(S, P) {
    P = fullParams(P); const g = P.geom;
    const halos = new Set(XB_HALOGENS); if (g.includeF) halos.add('F');
    const donors = [], accAtoms = [], accRings = [];
    const chainOK = a => !g.selectedChains || g.selectedChains.includes(a.chain);
    for (const a of S.atoms) {
      if (!chainOK(a)) continue; const R = S.residues[a.res];
      if (halos.has(a.el) && (R.cls === 'ligand' || (R.cls === 'protein' && R.het))) { const A = halogenPartner(S, a.i); if (A >= 0) donors.push({ X: a.i, A, res: a.res }); }
    }
    const scope = g.acceptorScope;
    for (const a of S.atoms) {
      if (!chainOK(a)) continue; const R = S.residues[a.res];
      const okCls = (R.cls === 'protein' && scope.protein) || (R.cls === 'nucleic' && scope.nucleic) || (R.cls === 'water' && scope.water) || (R.cls === 'ligand' && scope.otherLigand);
      if (!okCls) continue;
      if (['N', 'O', 'S', 'SE'].includes(a.el)) { const lw = lewisClass(S, a.i); if (g.excludeNonLewis && lw === 'weak') continue; accAtoms.push({ Y: a.i, res: a.res, lewis: lw }); }
    }
    if (g.includePi && scope.protein) for (const R of S.residues) {
      if (R.cls !== 'protein') continue; if (g.selectedChains && !g.selectedChains.includes(R.chain)) continue;
      ringsOf(S, R.idx, g.includeHisPi).forEach(rg => accRings.push({ res: R.idx, ids: rg.ids, label: rg.label }));
    }
    let maxCut = 0; Object.values(g.cutoffs).forEach(c => Object.values(c).forEach(v => { if (v > maxCut) maxCut = v; }));
    return { P, donors, accAtoms, accRings, maxCut };
  }
  function detectFrame(S, prep, frame, box) {
    const { P, donors, accAtoms, accRings, maxCut } = prep; const g = P.geom; const out = [];
    const mic = box && box[0] >= 10 && box[1] >= 10 && box[2] >= 10 && Math.abs(box[3] - 90) < 1e-3 && Math.abs(box[4] - 90) < 1e-3 && Math.abs(box[5] - 90) < 1e-3;
    const wrap = v => mic ? [v[0] - box[0] * Math.round(v[0] / box[0]), v[1] - box[1] * Math.round(v[1] / box[1]), v[2] - box[2] * Math.round(v[2] / box[2])] : v;
    const maxCut2 = maxCut * maxCut;
    for (const dn of donors) {
      const el = S.atoms[dn.X].el; const cut = g.cutoffs[el] || g.cutoffs.I;
      const px = coordOf(frame, dn.X), pa = coordOf(frame, dn.A);
      for (const ac of accAtoms) {
        if (ac.res === dn.res) continue;
        const vy = wrap(sub(coordOf(frame, ac.Y), px)); const d2 = dot(vy, vy); if (d2 > maxCut2) continue;
        const accEl = S.atoms[ac.Y].el; const c = cut[accEl === 'SE' ? 'S' : accEl]; if (!c) continue;
        const d = Math.sqrt(d2); if (d > c) continue;
        const py = add(px, vy); const th1 = angle(pa, px, py);
        const ok1 = th1 >= g.theta1Min && th1 <= g.theta1Max;
        let th2 = null; for (const r of S.bonds[ac.Y]) { if (S.atoms[r].el === 'H') continue; const t = angle(px, py, add(py, sub(coordOf(frame, r), coordOf(frame, ac.Y)))); th2 = th2 === null ? t : Math.min(th2, t); }
        const ok2 = !g.theta2Use || th2 === null || (th2 >= g.theta2Min && th2 <= g.theta2Max);
        if (!(ok1 && ok2) && !g.reportNearMiss) continue;
        const gs = geometricScore(el, accEl === 'SE' ? 'S' : accEl, d, th1, th2, P, false);
        out.push(Object.assign({ key: `${dn.X}>${ac.Y}`, type: 'atom', X: dn.X, A: dn.A, Y: ac.Y, donorRes: dn.res, accRes: ac.res, halogen: el, accEl, d, theta1: th1, theta2: th2,
          lewis: ac.lewis, status: ok1 && ok2 ? 'XB' : 'near-miss' }, gs, ok1 && ok2 ? {} : { G: 0 }));
      }
      for (const rg of accRings) {
        if (rg.res === dn.res) continue;
        const pts = rg.ids.map(i => add(px, wrap(sub(coordOf(frame, i), px))));
        const cen = centroid(pts); const v = sub(cen, px); const d = norm(v); const c = cut.AR; if (!c || d > c) continue;
        const th1 = angle(pa, px, cen); const nrm = planeNormal(pts);
        const piAng = deg(Math.acos(Math.min(1, Math.abs(dot(unit(v), nrm)))));
        const ok1 = th1 >= g.theta1Min && th1 <= g.theta1Max && piAng <= g.piNormalMax;
        if (!ok1 && !g.reportNearMiss) continue;
        const gs = geometricScore(el, 'AR', d, th1, null, P, true, piAng);
        out.push(Object.assign({ key: `${dn.X}>R${rg.ids[0]}`, type: 'pi', X: dn.X, A: dn.A, ring: rg.ids, ringLabel: rg.label, donorRes: dn.res, accRes: rg.res, halogen: el, accEl: 'π',
          d, theta1: th1, theta2: piAng, lewis: 'strong', status: ok1 ? 'XB' : 'near-miss' }, gs, ok1 ? {} : { G: 0 }));
      }
    }
    return out;
  }
  function detectXB(S, P, frameIdx = 0) {
    const prep = prepareDetection(S, P);
    const hits = detectFrame(S, prep, S.frames[frameIdx], S.boxes[frameIdx]);
    return { prep, hits: hits.map(h => decorate(S, h)).sort((a, b) => b.G - a.G || a.d - b.d) };
  }
  function decorate(S, h) {
    h.donorLabel = resLabel(S, h.donorRes); h.Xname = S.atoms[h.X].name; h.Aname = S.atoms[h.A].name;
    h.accLabel = resLabel(S, h.accRes); h.accName = h.type === 'pi' ? h.ringLabel : S.atoms[h.Y].name;
    h.accClass = h.type === 'pi' ? 'aromatic π' : acceptorClass(S, h.Y);
    return h;
  }

  /* ------------------------------------------------------------------ */
  /*  T2: trajectory persistence analysis                                */
  /* ------------------------------------------------------------------ */
  function parseXYZFrames(text, nAtoms) {
    const L = text.split(/\r?\n/); const frames = []; let i = 0;
    while (i < L.length) {
      const n = parseInt(L[i]); if (!n) { i++; continue; }
      if (nAtoms && n !== nAtoms) throw new Error(`XYZ frame has ${n} atoms but topology has ${nAtoms}`);
      if (i + 1 + n >= L.length + 1) throw new Error('XYZ file ends inside a frame');
      const f = new Float32Array(n * 3); const els = [];
      for (let k = 0; k < n; k++) { const t = (L[i + 2 + k] || '').trim().split(/\s+/); if (t.length < 4) throw new Error(`XYZ line ${i + 3 + k} is incomplete`); els.push(t[0].toUpperCase()); f[3 * k] = +t[1]; f[3 * k + 1] = +t[2]; f[3 * k + 2] = +t[3]; }
      frames.push({ f, els, comment: L[i + 1] }); i += n + 2;
    }
    return frames;
  }
  function parseDCD(buf, nAtomsExpected) {
    const dv = new DataView(buf); let le = true;
    let off = 0; let r = dv.getInt32(0, true); if (r !== 84) { le = false; r = dv.getInt32(0, false); if (r !== 84) throw new Error('Not a DCD file'); }
    const i32 = o => dv.getInt32(o, le);
    if (String.fromCharCode(dv.getUint8(4), dv.getUint8(5), dv.getUint8(6), dv.getUint8(7)) !== 'CORD') throw new Error('Not a DCD (CORD) file');
    const nset = i32(8); const nfixed = i32(4 + 4 + 8 * 4); const charmm = i32(4 + 80) !== 0; const hasUC = charmm && i32(4 + 4 + 10 * 4) === 1;
    off = 4 + 84 + 4; const tl = i32(off); off += 4 + tl + 4;
    off += 4; const natom = i32(off); off += 4 + 4;
    if (nfixed) throw new Error('DCD with fixed atoms not supported');
    if (nAtomsExpected && natom !== nAtomsExpected) throw new Error(`DCD has ${natom} atoms but topology has ${nAtomsExpected}`);
    const frames = []; const boxes = [];
    while (off < buf.byteLength - 4) {
      let box = null;
      if (hasUC) { const l = i32(off); off += 4; const u = []; for (let k = 0; k < 6; k++) u.push(dv.getFloat64(off + 8 * k, le)); off += l + 4; box = [u[0], u[2], u[5], 90, 90, 90]; }
      const f = new Float32Array(natom * 3);
      for (let c = 0; c < 3; c++) { const l = i32(off); off += 4; for (let k = 0; k < natom; k++) f[3 * k + c] = dv.getFloat32(off + 4 * k, le); off += l + 4; }
      frames.push(f); boxes.push(box);
      if (nset && frames.length >= nset && off >= buf.byteLength) break;
    }
    return { frames, boxes, natom };
  }
  function dropVsites(S, frames) {
    // frames read from a DCD/XYZ written with the virtual sites present -> frames in the atom order of S
    const v = S.vsites; if (!v) return frames;
    return frames.map(f => { if (f.length === S.atoms.length * 3) return f; if (f.length !== v.fileAtoms * 3) throw new Error(`frame has ${f.length / 3} atoms, topology ${S.atoms.length} (${v.fileAtoms} with virtual sites)`);
      const o = new Float32Array(v.keep.length * 3); v.keep.forEach((k, n) => { o[3 * n] = f[3 * k]; o[3 * n + 1] = f[3 * k + 1]; o[3 * n + 2] = f[3 * k + 2]; }); return o; });
  }
  function trajAccumulator(S, P, opts = {}) {
    P = fullParams(P); const prep = prepareDetection(S, P); const stride = opts.stride || P.dyn.stride || 1; const occCut = opts.occupancyCutoff ?? P.dyn.occCut;
    const stats = new Map(); let used = 0;
    return {
      add(frame, box, f) {
        used++;
        for (const h of detectFrame(S, prep, frame, box || null)) {
          if (h.status !== 'XB') continue;
          let s = stats.get(h.key); if (!s) { s = { proto: h, frames: [], d: [], t1: [], t2: [] }; stats.set(h.key, s); }
          s.frames.push(f); s.d.push(h.d); s.t1.push(h.theta1); s.t2.push(h.theta2);
        }
      },
      finish() {
        const res = [];
        for (const [key, s] of stats) {
          const occ = s.frames.length / used; if (occ * 100 < occCut) continue;
          let best = 0; for (let k = 1; k < s.d.length; k++) if (s.d[k] < s.d[best] - 1e-9 || (Math.abs(s.d[k] - s.d[best]) < 1e-9 && s.t1[k] > s.t1[best])) best = k;
          const runs = []; let run = 1;
          for (let k = 1; k < s.frames.length; k++) { if (s.frames[k] - s.frames[k - 1] === stride) run++; else { runs.push(run); run = 1; } }
          runs.push(run);
          const mean = a => a.reduce((x, y) => x + y, 0) / a.length; const sd = a => { const m = mean(a); return Math.sqrt(mean(a.map(v => (v - m) * (v - m)))); };
          const p = decorate(S, Object.assign({}, s.proto));
          const tauFrames = mean(runs) * stride; // in frames of the original trajectory
          res.push(Object.assign(p, { key, occupancy: occ, nFrames: s.frames.length, framesUsed: used, meanD: mean(s.d), sdD: sd(s.d), meanT1: mean(s.t1), sdT1: sd(s.t1),
            bestFrame: s.frames[best], bestD: s.d[best], bestT1: s.t1[best], meanResidence: tauFrames, maxResidence: Math.max(...runs) * stride, nEvents: runs.length,
            series: { frames: s.frames, d: s.d, t1: s.t1 }, D: Dscore(occ, tauFrames, P) }));
        }
        res.sort((a, b) => b.occupancy - a.occupancy);
        return { results: res, framesUsed: used };
      }
    };
  }
  function Dscore(occ, tauFrames, P) { const d = (P || DEFAULTS).dyn; return occ * Math.pow(1 - Math.exp(-tauFrames / d.tauScale), d.dExp); }
  function analyzeTrajectory(S, frames, boxes, P, opts = {}) {
    P = fullParams(P); const acc = trajAccumulator(S, P, opts); const stride = opts.stride || P.dyn.stride || 1;
    for (let f = 0; f < frames.length; f += stride) acc.add(frames[f], boxes ? boxes[f] : null, f);
    const r = acc.finish(); r.totalFrames = frames.length; return r;
  }

  /* ------------------------------------------------------------------ */
  /*  Materials: CIF reader with symmetry expansion + periodic contacts  */
  /* ------------------------------------------------------------------ */
  function parseCIF(text) {
    const L = text.split(/\r?\n/); const num = s => parseFloat(String(s).replace(/\(.*\)/, ''));
    const cell = {}; const ops = []; const loops = [];
    for (let i = 0; i < L.length; i++) {
      const t = L[i].trim(); let m;
      if ((m = /^_cell_length_([abc])\s+(\S+)/.exec(t))) cell[m[1]] = num(m[2]);
      if ((m = /^_cell_angle_(alpha|beta|gamma)\s+(\S+)/.exec(t))) cell[m[1]] = num(m[2]);
      if (t === 'loop_') {
        const heads = []; let j = i + 1;
        while (j < L.length && L[j].trim().startsWith('_')) { heads.push(L[j].trim().split(/\s+/)[0]); j++; }
        const rows = [];
        while (j < L.length) { const s = L[j].trim(); if (!s || s.startsWith('_') || s === 'loop_' || s.startsWith('data_') || s.startsWith('#')) break; rows.push(s.match(/'[^']*'|"[^"]*"|\S+/g).map(x => x.replace(/^['"]|['"]$/g, ''))); j++; }
        loops.push({ heads, rows }); i = j - 1;
      }
    }
    for (const lp of loops) { const k = lp.heads.findIndex(h => /_symmetry_equiv_pos_as_xyz|_space_group_symop_operation_xyz/.test(h)); if (k >= 0) lp.rows.forEach(r => ops.push(r.length > lp.heads.length ? r.slice(k).join('') : r[k])); }
    if (!ops.length) ops.push('x,y,z');
    const site = loops.find(lp => lp.heads.includes('_atom_site_fract_x')); if (!site) throw new Error('CIF has no _atom_site_fract_x loop');
    ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].forEach(k => { if (!(cell[k] > 0)) throw new Error('CIF cell parameter _cell_' + k + ' missing'); });
    const H = n => site.heads.indexOf(n);
    const iL = H('_atom_site_label'), iT = H('_atom_site_type_symbol'), ix = H('_atom_site_fract_x'), iy = H('_atom_site_fract_y'), iz = H('_atom_site_fract_z');
    const asym = site.rows.filter(r => r.length >= site.heads.length).map(r => ({ label: r[iL], el: (iT >= 0 ? r[iT] : r[iL]).replace(/[^A-Za-z].*$/, '').toUpperCase(), f: [num(r[ix]), num(r[iy]), num(r[iz])] }));
    const fns = ops.map(op => op.toLowerCase().replace(/\s/g, '').split(',').map(e => { if (!/^[xyz0-9+\-*/.()]+$/.test(e)) throw new Error('Unsupported symmetry operator: ' + op); return new Function('x', 'y', 'z', 'return ' + e.replace(/(\d)\/(\d)/g, '($1/$2)') + ';'); }));
    const frac = [];
    for (const a of asym) for (const fn of fns) {
      const p = fn.map(g => g(a.f[0], a.f[1], a.f[2])).map(v => v - Math.floor(v));
      if (!frac.some(q => q.el === a.el && q.f.every((v, k) => { const dd = Math.abs(v - p[k]); return Math.min(dd, 1 - dd) < 1e-3; }))) frac.push({ label: a.label, el: a.el, f: p });
    }
    const M = cellMatrix(cell);
    return { cell, M, ops, asymCount: asym.length, atoms: frac.map(a => Object.assign(a, { c: fracToCart(M, a.f) })) };
  }
  function cellMatrix(c) { const al = c.alpha * Math.PI / 180, be = c.beta * Math.PI / 180, ga = c.gamma * Math.PI / 180;
    const ax = [c.a, 0, 0], bx = [c.b * Math.cos(ga), c.b * Math.sin(ga), 0]; const cx = c.c * Math.cos(be), cy = c.c * (Math.cos(al) - Math.cos(be) * Math.cos(ga)) / Math.sin(ga);
    return [ax, bx, [cx, cy, Math.sqrt(Math.max(0, c.c * c.c - cx * cx - cy * cy))]]; }
  function fracToCart(M, f) { return [f[0] * M[0][0] + f[1] * M[1][0] + f[2] * M[2][0], f[0] * M[0][1] + f[1] * M[1][1] + f[2] * M[2][1], f[0] * M[0][2] + f[1] * M[1][2] + f[2] * M[2][2]]; }
  function parseXYZ(text) { const fr = parseXYZFrames(text)[0]; if (!fr) throw new Error('No XYZ frame found'); return { atoms: fr.els.map((el, k) => ({ label: el + (k + 1), el, c: [fr.f[3 * k], fr.f[3 * k + 1], fr.f[3 * k + 2]] })) }; }
  function imageList(model, n) { const imgs = []; if (model.M) { for (let i = -n; i <= n; i++) for (let j = -n; j <= n; j++) for (let k = -n; k <= n; k++) imgs.push([i, j, k]); } else imgs.push([0, 0, 0]); return imgs; }
  function detectMaterialXB(model, P) {
    P = fullParams(P); const g = P.geom;
    const at = model.atoms; const imgs = imageList(model, P.mat.images);
    const shift = t => model.M ? fracToCart(model.M, t) : [0, 0, 0];
    const out = []; const halos = new Set(XB_HALOGENS); if (g.includeF) halos.add('F');
    at.forEach((x, xi) => {
      if (!halos.has(x.el)) return;
      let A = null, Ad = 9;
      at.forEach((b, bi) => imgs.forEach(t => { if (bi === xi && t.every(v => v === 0)) return; const p = add(b.c, shift(t)); const d = dist(x.c, p); if (d < Ad && d < (COV[x.el] + (COV[b.el] || 0.8) + P.struct.bondTol)) { Ad = d; A = { bi, p }; } }));
      if (!A) return;
      at.forEach((y, yi) => {
        if (!['N', 'O', 'S'].includes(y.el)) return;
        imgs.forEach(t => {
          const p = add(y.c, shift(t)); const d = dist(x.c, p); const c = (g.cutoffs[x.el] || {})[y.el];
          if (!c || d > c || d < P.mat.minDist) return;
          const th1 = angle(A.p, x.c, p); const ok = th1 >= g.theta1Min && th1 <= g.theta1Max;
          if (!ok && !g.reportNearMiss) return;
          const gs = geometricScore(x.el, y.el, d, th1, null, P, false);
          out.push(Object.assign({ X: xi, Xlabel: x.label, A: A.bi, Alabel: at[A.bi].label, Y: yi, Ylabel: y.label, image: t.join(' '), halogen: x.el, accEl: y.el, d, theta1: th1, status: ok ? 'XB' : 'near-miss' }, gs, ok ? {} : { G: 0 }, { Ypos: p, Apos: A.p }));
        });
      });
    });
    // Symmetry copies of one contact share their atom labels and geometry. They are reported once with a multiplicity.
    const uniq = new Map();
    out.forEach(h => { const k = [h.Alabel, h.Xlabel, h.Ylabel, h.d.toFixed(3), h.theta1.toFixed(1), h.status].join('|');
      if (uniq.has(k)) uniq.get(k).mult++; else { h.mult = 1; uniq.set(k, h); } });
    return [...uniq.values()].sort((a, b) => b.G - a.G);
  }
  function buildMaterialQM(model, hit, opts = {}) {
    const P = fullParams(opts.P); const at = model.atoms; const imgs = imageList(model, P.mat.images);
    const shift = t => model.M ? fracToCart(model.M, t) : [0, 0, 0];
    function grow(startIdx, startPos) {
      const mol = [{ i: startIdx, p: startPos }]; const seen = new Map([[startIdx, startPos]]);
      for (let q = 0; q < mol.length && mol.length < 400; q++) {
        const cur = mol[q];
        at.forEach((b, bi) => { if (seen.has(bi)) return;
          for (const t of imgs) { const pb = add(b.c, shift(t)); const d = dist(cur.p, pb);
            if ((cur.i === hit.X && bi === hit.Y) || (cur.i === hit.Y && bi === hit.X)) continue; // the halogen bond itself is never a covalent bond
            if (d > 0.4 && d < (COV[at[cur.i].el] || 0.8) + (COV[b.el] || 0.8) + P.struct.bondTol) { seen.set(bi, pb); mol.push({ i: bi, p: pb }); break; } } });
      }
      return mol;
    }
    const m1 = grow(hit.X, at[hit.X].c); const m2 = grow(hit.Y, hit.Ypos);
    const same = m1.some(a => a.i === hit.Y && dist(a.p, hit.Ypos) < 1e-3);
    const qm = []; const indexOf1 = new Map(), indexOf2 = new Map();
    m1.forEach(a => { indexOf1.set(a.i, qm.length); qm.push({ el: at[a.i].el, x: a.p[0], y: a.p[1], z: a.p[2], frag: 1, src: { label: at[a.i].label } }); });
    if (!same) m2.forEach(a => { indexOf2.set(a.i, qm.length); qm.push({ el: at[a.i].el, x: a.p[0], y: a.p[1], z: a.p[2], frag: 2, src: { label: at[a.i].label + "'" } }); });
    const warnings = [];
    if (!qm.some(a => a.el === 'H')) warnings.push('No hydrogens in the structure: add H (e.g. riding model) before QM; H positions are not generated for periodic inputs.');
    if (same) warnings.push('Halogen and acceptor belong to the same molecule instance (intramolecular contact); EDA fragments cannot be defined.');
    if (m1.length >= 400 || m2.length >= 400) warnings.push('Fragment growth hit 400 atoms: polymeric / framework material. Truncate manually.');
    const fragCharge = Object.assign({ 1: 0, 2: 0 }, opts.chargeOverride || {});
    const totalCharge = fragCharge[1] + fragCharge[2];
    const electrons = qm.reduce((s, a) => s + (Z[a.el] || 0), 0) - totalCharge;
    if (electrons % 2) warnings.push(`Odd electron count (${electrons}); set fragment charges/multiplicity.`);
    const xbs = [{ id: hit.cid || 'XB1', type: 'atom', X: indexOf1.get(hit.X) + 1, A: indexOf1.get(hit.A) + 1, Y: (same ? indexOf1.get(hit.Y) : indexOf2.get(hit.Y)) + 1, ring: null,
      label: `${hit.Xlabel} ··· ${hit.Ylabel}' [${hit.image}]`, d: hit.d, theta1: hit.theta1, halogen: hit.halogen }];
    const manifest = { tool: 'BondX', version: VERSION, created: new Date().toISOString(), source: opts.sourceName || 'material', frame: 0, mode: 'crystal-dimer', shell: 0,
      nAtoms: qm.length, addedH: 0, linkH: 0, fragCharge, totalCharge, multiplicity: opts.multiplicity || (electrons % 2 ? 2 : 1), xbs,
      atoms: qm.map((a, k) => ({ n: k + 1, el: a.el, frag: a.frag, x: +a.x.toFixed(5), y: +a.y.toFixed(5), z: +a.z.toFixed(5), src: a.src.label })), warnings,
      assumptions: ['Crystal dimer: the two molecules forming the contact are grown from covalent connectivity across periodic images.', 'Both fragments are neutral unless overridden.'] };
    manifest.id = 'BX-' + hashString(JSON.stringify(manifest.atoms) + JSON.stringify(xbs)).slice(0, 10);
    return { qm, manifest, include: new Set(), indexOf: new Map() };
  }

  /* ------------------------------------------------------------------ */
  /*  T3: force-field typing (protonation states + partial charges)      */
  /* ------------------------------------------------------------------ */
  const HEAVY_ALIAS = { OT1: 'O', OT2: 'OXT', O1P: 'OP1', O2P: 'OP2', C5M: 'C7', OW: 'O', OH2: 'O' };
  const H_ALIAS = { HN: 'H', HT1: 'H1', HT2: 'H2', HT3: 'H3', HN1: 'H1', HN2: 'H2', HW1: 'H1', HW2: 'H2' };
  const CHARMM_ION = { NA: 'SOD', K: 'POT', CL: 'CLA', CA: 'CAL', MG: 'MG', ZN: 'ZN2', CD: 'CD2', LI: 'LIT', RB: 'RUB', CS: 'CES', BA: 'BAR' };
  const ION_FORMAL = { NA: 1, K: 1, LI: 1, RB: 1, CS: 1, TL: 1, AG: 1, CU1: 1, MG: 2, CA: 2, ZN: 2, MN: 2, FE2: 2, CO: 2, NI: 2, CU: 2, CD: 2, HG: 2, SR: 2, BA: 2, PB: 2, PD: 2, PT: 2,
    FE: 3, AL: 3, CR: 3, GA: 3, IN: 3, LA: 3, EU3: 3, GD3: 3, CL: -1, BR: -1, IOD: -1, I: -1, F: -1 };
  function ffSet(name) { const f = FFDATA && FFDATA[name]; if (!f) throw new Error('Unknown force field ' + name); return f; }
  function resId(R) { return `${R.chain}:${R.resSeq}${R.iCode || ''}`; }
  // Choose the force-field template name for residue R
  function ffVariant(S, R, P) {
    const f = P.ff; const FF = ffSet(f.name); const ov = (f.overrides || {})[resId(R)];
    if (R.cls === 'water') return { tpl: 'WATER', why: 'water model ' + f.water };
    if (R.cls === 'ion') return { tpl: 'ION', why: 'ion' };
    if (R.cls === 'nucleic') {
      let base = R.resName; if (/^R[ACGU]$/.test(base)) base = base.slice(1);
      const five = R.nterm && R.byName.P === undefined, three = !!R.cterm;
      const key = base + (five && three ? 'N' : five ? '5' : three ? '3' : '');
      return FF.residues[key] ? { tpl: key, why: five ? "5'-terminal" : three ? "3'-terminal" : 'internal' } : { tpl: null, why: 'no template ' + key };
    }
    if (R.cls !== 'protein') return { tpl: null, why: 'non-standard residue' };
    const has = n => R.byName[n] !== undefined; let base, why;
    const rn = R.resName;
    if (ov) { base = ov; why = 'user choice'; }
    else if (rn === 'HIS' || rn === 'HID' || rn === 'HIE' || rn === 'HIP' || rn === 'HSD' || rn === 'HSE' || rn === 'HSP') {
      const named = { HID: 'HID', HIE: 'HIE', HIP: 'HIP', HSD: 'HID', HSE: 'HIE', HSP: 'HIP' }[rn];
      if (named) { base = named; why = 'named in file'; }
      else if (has('HD1') && has('HE2')) { base = 'HIP'; why = 'HD1 and HE2 present'; }
      else if (has('HD1')) { base = 'HID'; why = 'HD1 present'; } else if (has('HE2')) { base = 'HIE'; why = 'HE2 present'; }
      else { base = f.hisDefault; why = 'default'; }
    } else if (rn === 'CYS' || rn === 'CYX' || rn === 'CYM') {
      if (rn !== 'CYS') { base = rn; why = 'named in file'; }
      else {
        const sg = R.byName.SG; let ss = false;
        if (sg !== undefined) { const p = coordOf(S.frames[0], sg); for (const R2 of S.residues) { if (R2 === R || R2.cls !== 'protein' || AA_ALIAS[R2.resName] !== 'CYS' && R2.resName !== 'CYS') continue; const s2 = R2.byName.SG; if (s2 !== undefined && dist(p, coordOf(S.frames[0], s2)) < f.ssCut) { ss = true; break; } } }
        base = ss ? 'CYX' : 'CYS'; why = ss ? `S–S < ${f.ssCut} Å` : 'free thiol';
      }
    } else if (rn === 'ASP' || rn === 'ASH') { base = rn === 'ASH' || has('HD2') ? 'ASH' : 'ASP'; why = base === 'ASH' ? 'protonated in file' : 'pH 7 rule'; }
    else if (rn === 'GLU' || rn === 'GLH') { base = rn === 'GLH' || has('HE2') ? 'GLH' : 'GLU'; why = base === 'GLH' ? 'protonated in file' : 'pH 7 rule'; }
    else if (rn === 'LYS' || rn === 'LYN') { base = rn === 'LYN' || (has('HZ1') && has('HZ2') && !has('HZ3')) ? 'LYN' : 'LYS'; why = base === 'LYN' ? 'neutral in file' : 'pH 7 rule'; }
    else if (rn === 'MSE') { base = 'MET'; why = 'selenomethionine treated as MET (Se takes the SD charge)'; }
    else { base = rn; why = 'standard'; }
    const isN = R.nterm === 'start' || (R.nterm === 'break' && f.chainBreaks === 'charged');
    const isC = R.cterm === 'end' || (R.cterm === 'break' && f.chainBreaks === 'charged');
    let tpl = base;
    if (f.termini === 'charged' && (isN || isC)) {
      const pre = isN ? 'N' : 'C';
      if (FF.residues[pre + base]) { tpl = pre + base; why += isN ? '; N-terminus' : '; C-terminus'; if (isN && isC) why += ' (single-residue chain: C-terminal charge not added)'; }
      else why += `; no ${pre}-terminal template for ${base} (internal charges used)`;
    }
    if (!FF.residues[tpl]) return { tpl: null, why: 'no template ' + tpl };
    return { tpl, why, base };
  }
  // Assign a partial charge to every atom (and a list of folded hydrogens) — the single source for QM and MM
  function assignForceField(S, P, extra = {}) {
    P = fullParams(P); const f = P.ff; const FF = ffSet(f.name); const frame = S.frames[0];
    const q = new Float64Array(S.atoms.length); const type = new Array(S.atoms.length).fill('');
    const res = []; const warnings = []; let usePQR = f.chargeSource === 'pqr' || (f.chargeSource === 'auto' && S.isPQR);
    if (f.chargeSource === 'pqr' && !S.isPQR) { warnings.push('PQR charges requested but the structure is not a PQR file; force-field charges used.'); usePQR = false; }
    const mol2 = extra.mol2 || null; // {resName: {atomName: q}}  (hydrogens already folded where needed)
    const waterT = FF.waters[f.water] || Object.values(FF.waters)[0];
    for (const R of S.residues) {
      const info = { idx: R.idx, id: resId(R), resName: R.resName, cls: R.cls, tpl: null, why: '', net: 0, formal: null, missing: [], extra: [], foldedH: 0, source: '' };
      if (usePQR) {
        R.atoms.forEach(i => { q[i] = S.atoms[i].q || 0; }); info.source = 'PQR'; info.net = R.atoms.reduce((s, i) => s + q[i], 0); info.formal = Math.round(info.net);
        const v = R.cls === 'protein' ? ffVariant(S, R, P) : { tpl: null }; info.tpl = v.tpl; info.why = 'charges from PQR file';
        res.push(info); continue;
      }
      if (R.cls === 'ion') {
        const a = S.atoms[R.atoms[0]]; const key = f.name === 'charmm36' ? (CHARMM_ION[R.resName] || R.resName) : R.resName;
        let c = FF.ions[key]; let src = f.name;
        if (c === undefined) { c = ION_FORMAL[R.resName] ?? ION_FORMAL[a.el]; src = 'formal'; }
        if (c === undefined) { c = 0; warnings.push(`Ion ${R.resName} ${info.id}: unknown charge, set to 0.`); src = 'unknown'; }
        q[a.i] = c; info.net = c; info.formal = Math.round(c); info.tpl = key; info.why = src === 'formal' ? 'formal charge' : 'ion parameters'; info.source = src; res.push(info); continue;
      }
      let tpl = null;
      if (R.cls === 'water') { tpl = waterT; info.tpl = f.water; info.why = 'water model'; }
      else if (R.cls === 'protein' || R.cls === 'nucleic') { const v = ffVariant(S, R, P); info.tpl = v.tpl; info.why = v.why; tpl = v.tpl ? FF.residues[v.tpl] : null; }
      if (!tpl) {
        if (mol2 && mol2[R.resName]) {
          const m = mol2[R.resName]; let missing = 0;
          R.atoms.forEach(i => { const a = S.atoms[i]; if (m[a.name] !== undefined) q[i] = m[a.name]; else { q[i] = 0; if (a.el !== 'H') missing++; } });
          info.source = 'MOL2'; info.net = R.atoms.reduce((s, i) => s + q[i], 0); info.formal = Math.round(m._total ?? info.net); info.why = 'charges from MOL2';
          if (missing) warnings.push(`${R.resName} ${info.id}: ${missing} atom(s) not found in the MOL2 file (charge 0).`);
        } else {
          // no template and no MOL2: the formal charge from group rules is spread evenly over the heavy atoms
          // (e.g. sulfate -2 -> -0.4 on S and each O); neutral groups stay at 0
          const lpq = R.atoms.some(i => S.atoms[i].el === 'H') ? null : ligandProtonation(S, S.frames[0], R, 'rules');
          const g = lpq && !lpq.metal ? lpq.charge : ligandChargeGuess(S, S.frames[0], R); const heavy = R.atoms.filter(i => S.atoms[i].el !== 'H');
          if (g !== 0 && heavy.length) { R.atoms.forEach(i => { q[i] = S.atoms[i].el === 'H' ? 0 : g / heavy.length; }); info.source = 'formal'; info.net = g; info.formal = g; info.why = `no template or MOL2: formal charge ${g > 0 ? '+' : ''}${g} (group rules) spread over ${heavy.length} heavy atoms`; }
          else { info.source = 'none'; info.why = info.why || 'no force-field template; charges set to 0'; info.formal = g === 0 && heavy.length ? 0 : null; }
        }
        res.push(info); continue;
      }
      // template matching (heavy atoms by name, hydrogens by parent)
      info.source = f.name; const T = tpl.atoms; const used = new Set();
      const heavy = R.atoms.filter(i => S.atoms[i].el !== 'H');
      const nameOf = a => { let n = a.name; if (R.resName === 'MSE' && n === 'SE') n = 'SD'; if (T[n] === undefined && HEAVY_ALIAS[n] && T[HEAVY_ALIAS[n]] !== undefined) n = HEAVY_ALIAS[n];
        if (T[n] === undefined && n.includes('*')) n = n.replace(/\*/g, "'"); if (T[n] === undefined && R.resName === 'ILE' && n === 'CD') n = 'CD1'; return n; };
      for (const i of heavy) {
        const a = S.atoms[i]; const n = nameOf(a); const t = T[n];
        if (!t || used.has(n)) { q[i] = 0; info.extra.push(a.name); continue; }
        used.add(n); q[i] = t[0]; type[i] = t[1];
        // hydrogens bonded to this heavy atom
        const hs = S.bonds[i].filter(j => S.atoms[j].el === 'H' && S.atoms[j].res === R.idx);
        const tH = t[2].slice(); const assigned = new Set();
        for (const j of hs) { const hn = S.atoms[j].name; let k = tH.findIndex(([nm], kk) => !assigned.has(kk) && (nm === hn || nm === H_ALIAS[hn] || H_ALIAS[nm] === hn));
          if (k < 0) k = tH.findIndex((_, kk) => !assigned.has(kk));
          if (k < 0) { q[j] = 0; info.extra.push(hn); continue; }
          assigned.add(k); q[j] = tH[k][1]; type[j] = 'H'; }
        tH.forEach(([nm, qh], kk) => { if (assigned.has(kk)) return; if (f.foldH) { q[i] += qh; info.foldedH++; } });
      }
      R.atoms.filter(i => S.atoms[i].el === 'H' && !S.bonds[i].some(j => S.atoms[j].res === R.idx && S.atoms[j].el !== 'H')).forEach(j => { q[j] = 0; info.extra.push(S.atoms[j].name); });
      // missing heavy atoms
      let missQ = 0; for (const [n, t] of Object.entries(T)) if (!used.has(n)) { info.missing.push(n); missQ += t[0] + (f.foldH ? t[2].reduce((s, h) => s + h[1], 0) : 0); }
      if (info.missing.length && f.missing === 'redistribute') { const present = heavy.filter(i => type[i]); if (present.length) present.forEach(i => { q[i] += missQ / present.length; }); }
      info.net = R.atoms.reduce((s, i) => s + q[i], 0); info.formal = R.cls === 'nucleic' ? Math.round(tpl.q) : tpl.q;
      if (info.missing.length) warnings.push(`${R.resName} ${info.id}: ${info.missing.length} atom(s) missing (${info.missing.slice(0, 6).join(', ')}${info.missing.length > 6 ? '…' : ''})` + (f.missing === 'redistribute' ? `; their charge (${missQ.toFixed(3)} e) spread over the residue.` : '; net charge changed.'));
      if (info.extra.length) warnings.push(`${R.resName} ${info.id}: atom(s) not in the ${info.tpl} template get charge 0 (${info.extra.slice(0, 6).join(', ')}).`);
      res.push(info);
    }
    const total = q.reduce((s, v) => s + v, 0);
    return { q, type, residues: res, warnings, total, ff: f.name, water: f.water, usedPQR: usePQR };
  }
  // titratable residues and their current variants (for the protonation table)
  function titratableResidues(S, P, ffa) {
    P = fullParams(P); const out = [];
    S.residues.forEach(R => { if (R.cls !== 'protein') return; const b = AA_ALIAS[R.resName] || R.resName;
      if (!['HIS', 'ASP', 'GLU', 'LYS', 'CYS', 'TYR', 'ARG'].includes(b)) return; const v = ffVariant(S, R, P); const info = ffa ? ffa.residues[R.idx] : null;
      out.push({ idx: R.idx, id: resId(R), resName: R.resName, base: b, variant: v.base || v.tpl, tpl: v.tpl, why: v.why, net: info ? info.net : null,
        choices: { HIS: ['HIE', 'HID', 'HIP'], ASP: ['ASP', 'ASH'], GLU: ['GLU', 'GLH'], LYS: ['LYS', 'LYN'], CYS: ['CYS', 'CYX', 'CYM'], TYR: ['TYR'], ARG: ['ARG'] }[b] }); });
    return out;
  }
  // MOL2 parser: per residue name, atom-name -> charge, hydrogens folded onto heavy atoms when foldH
  function parseMOL2(text, foldH = true) {
    const sec = name => { const m = new RegExp('@<TRIPOS>' + name + '\\s*\\n([\\s\\S]*?)(?=@<TRIPOS>|$)').exec(text); return m ? m[1].trim().split(/\r?\n/).filter(l => l.trim()) : []; };
    const atoms = sec('ATOM').map(l => { const t = l.trim().split(/\s+/); return { id: +t[0], name: t[1], x: +t[2], y: +t[3], z: +t[4], type: t[5], resName: (t[7] || 'LIG').replace(/\d+$/, '').toUpperCase().slice(0, 4), q: +t[8], el: t[5].split('.')[0].toUpperCase() }; });
    if (!atoms.length || atoms.some(a => isNaN(a.q))) throw new Error('MOL2 file has no @<TRIPOS>ATOM section with charges');
    const mol = sec('MOLECULE'); if (mol[3] && /NO_CHARGES/i.test(mol[3])) throw new Error('MOL2 file declares NO_CHARGES: generate charges first (e.g. antechamber -c bcc for AM1-BCC/GAFF2).');
    if (atoms.every(a => a.q === 0)) throw new Error('All MOL2 charges are zero: generate charges first (e.g. antechamber -c bcc).');
    const bonds = sec('BOND').map(l => l.trim().split(/\s+/).slice(1, 3).map(Number));
    const byId = new Map(atoms.map(a => [a.id, a]));
    const parent = new Map(); bonds.forEach(([a, b]) => { const A = byId.get(a), Bb = byId.get(b); if (!A || !Bb) return; if (A.el === 'H' && Bb.el !== 'H') parent.set(A.id, Bb); if (Bb.el === 'H' && A.el !== 'H') parent.set(Bb.id, A); });
    const out = {};
    for (const a of atoms) { const r = out[a.resName] || (out[a.resName] = { _total: 0, _H: {} }); r._total += a.q; if (a.el !== 'H') r[a.name] = (r[a.name] || 0) + a.q; else r._H[a.name] = a.q; }
    for (const a of atoms) if (a.el === 'H') { const p = parent.get(a.id); const r = out[a.resName]; r[a.name] = a.q; if (foldH && p) r['__fold_' + p.name] = (r['__fold_' + p.name] || 0) + a.q; }
    // final map used by the engine: heavy atoms carry folded H charge (H present in a structure get their own charge)
    const res = {};
    for (const [rn, r] of Object.entries(out)) { res[rn] = { _total: Math.round(r._total * 1e4) / 1e4, _nAtoms: atoms.filter(a => a.resName === rn).length };
      for (const [k, v] of Object.entries(r)) { if (k.startsWith('_')) continue; res[rn][k] = v; }
      for (const [k, v] of Object.entries(r)) if (k.startsWith('__fold_')) res[rn]['__H_' + k.slice(7)] = v; }
    return res;
  }
  // structure has H? -> use per-atom; else fold. Applies the folded variant for a structure without H.
  function mol2ForStructure(mol2, S) {
    if (!mol2) return null; const out = {};
    for (const [rn, m] of Object.entries(mol2)) { out[rn] = { _total: m._total };
      const hasH = S.residues.some(R => R.resName === rn && R.atoms.some(i => S.atoms[i].el === 'H'));
      for (const [k, v] of Object.entries(m)) { if (k.startsWith('_')) continue; out[rn][k] = v + (!hasH && m['__H_' + k] ? m['__H_' + k] : 0); } }
    return out;
  }

  /* ------------------------------------------------------------------ */
  /*  T3a: hydrogen saturation helpers                                   */
  /* ------------------------------------------------------------------ */
  const SP2_ATOMS = { _BB: ['C', 'O', 'N'], ARG: ['NE', 'CZ', 'NH1', 'NH2'], ASN: ['CG', 'OD1', 'ND2'], ASP: ['CG', 'OD1', 'OD2'], GLN: ['CD', 'OE1', 'NE2'],
    GLU: ['CD', 'OE1', 'OE2'], HIS: ['CG', 'ND1', 'CD2', 'CE1', 'NE2'], PHE: ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    TRP: ['CG', 'CD1', 'NE1', 'CD2', 'CE2', 'CE3', 'CZ3', 'CZ2', 'CH2'], TYR: ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'] };
  function placeH(center, nbrPos, nH, hyb, ref, len) {
    const c = center; const u = nbrPos.map(p => unit(sub(p, c))); const out = []; const mk = d => add(c, scl(unit(d), len)); const n = u.length;
    if (n === 0) { const d0 = [1, 0, 0]; for (let k = 0; k < nH; k++) out.push(mk(k === 0 ? d0 : rotAround(d0, [0, 0, 1], k * (hyb === 'sp3' ? 109.47 : 120) * Math.PI / 180))); return out; }
    if (hyb === 'sp' && n === 1) { out.push(mk(scl(u[0], -1))); return out; }
    if (n >= 3) { out.push(mk(scl(add(add(u[0], u[1]), u[2]), -1))); return out; }
    if (n === 2) {
      const b = unit(scl(add(u[0], u[1]), -1));
      if (hyb === 'sp2' || (nH === 1 && hyb !== 'sp3')) { out.push(mk(b)); return out; }
      const p = unit(cross(u[0], u[1])); const h = 54.75 * Math.PI / 180;
      out.push(mk(add(scl(b, Math.cos(h)), scl(p, Math.sin(h))))); if (nH > 1) out.push(mk(add(scl(b, Math.cos(h)), scl(p, -Math.sin(h)))));
      return out;
    }
    const ax = scl(u[0], -1); let perp;
    if (ref) { const r = sub(ref, nbrPos[0]); perp = unit(sub(r, scl(u[0], dot(r, u[0])))); if (norm(perp) < 1e-6) perp = anyPerp(ax); } else perp = anyPerp(ax);
    if (hyb === 'sp2') { const s = Math.sin(60 * Math.PI / 180), co = Math.cos(60 * Math.PI / 180); out.push(mk(add(scl(ax, co), scl(perp, -s)))); if (nH > 1) out.push(mk(add(scl(ax, co), scl(perp, s)))); return out; }
    const tilt = 70.5 * Math.PI / 180;
    for (let k = 0; k < nH; k++) { const d0 = add(scl(ax, Math.cos(tilt)), scl(perp, -Math.sin(tilt))); out.push(mk(rotAround(d0, ax, k * 2 * Math.PI / 3))); }
    return out;
  }
  function guessHyb(el, center, nbrPos, nbrEls) {
    const n = nbrPos.length;
    if (n >= 2) {
      let sum = 0, cnt = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { sum += angle(nbrPos[i], center, nbrPos[j]); cnt++; }
      const avg = sum / cnt;
      if (n === 2 && avg > 165) return 'sp'; if (n === 3) return sum >= 350 ? 'sp2' : 'sp3';
      if (n === 2) { if (el === 'C' || el === 'N') { const dmin = Math.min(...nbrPos.map(p => dist(p, center))); return (avg > 114 || dmin < 1.40) ? 'sp2' : 'sp3'; } return 'sp3'; }
      return 'sp3';
    }
    if (n === 1) {
      const d = dist(nbrPos[0], center); const o = nbrEls[0];
      if (el === 'O') return (o === 'C' && d < 1.30) || (o === 'S' && d < 1.50) || (o === 'P' && d < 1.55) || (o === 'N' && d < 1.28) ? 'sp2' : 'sp3';
      if (el === 'N') return d < 1.20 ? 'sp' : (d < 1.39 ? 'sp2' : 'sp3');
      if (el === 'C') return d < 1.22 ? 'sp' : (d < 1.42 ? 'sp2' : 'sp3');
      if (el === 'S') return (o === 'C' && d < 1.70) || o === 'P' ? 'sp2' : 'sp3';
    }
    return 'sp3';
  }
  function heuristicHCount(el, hyb, nHeavy, center, nbrPos) {
    if (el === 'C') return Math.max(0, (hyb === 'sp3' ? 4 : hyb === 'sp2' ? 3 : 2) - nHeavy);
    if (el === 'N') { if (hyb === 'sp' || nHeavy >= 3) return 0; if (hyb === 'sp2' && nHeavy === 2) { const ds = nbrPos.map(p => dist(p, center)); return ds.every(d => d < 1.36) ? 0 : 1; } return Math.max(0, 3 - nHeavy); }
    if (el === 'O') { if (hyb === 'sp2') return 0; return Math.max(0, 2 - nHeavy); }
    if (el === 'S') { if (nHeavy >= 2 || hyb === 'sp2') return 0; return 1; }
    if (el === 'SE') return nHeavy >= 2 ? 0 : 1;
    return 0;
  }
  // ---------------------------------------------------------------------------------------------
  // Ligand protonation from heavy atoms (one consistent structure for H count and charge).
  // 1. sp2 candidates: C/N in planar unsaturated rings, planar 3-coordinate C, 2-coordinate C/N with
  //    sp2 geometry, terminal C=O/C=S/C=N by bond length. 3-coordinate N and 2-coordinate O/S never
  //    take a double bond (amide/pyrrole N, furan O).
  // 2. Kekule structure: every sp2 C gets exactly one double bond; among such structures the one with
  //    the most C=O/C=S (keto form) and then the shortest double bonds is chosen.
  // 3. H count = valence - heavy neighbours - double bond (C 4, N 3, O/S 2): an sp2 N left without a
  //    double bond is pyrrole-type and gets one H.
  // 4. P keeps one =O, S keeps up to two =O (the shortest bonds); other terminal O are acid OH.
  // 5. pH-7 group rules (mode 'rules'): COOH -> COO-, phosphate/phosphonate/phosphinate -1 or -2,
  //    sulfonate -1, sulfate -2 (the matching O-H are removed); R4N+ +1. Mode 'zero': neutral.
  function ligandProtonation(S, frame, R, mode) {
    const heavy = R.atoms.filter(i => S.atoms[i].el !== 'H'); const inR = new Set(heavy);
    const EL = i => S.atoms[i].el; const pos = i => coordOf(frame, i);
    const metal = heavy.find(i => !['C', 'N', 'O', 'S', 'P', 'F', 'CL', 'BR', 'I', 'B', 'SE', 'SI', 'AS'].includes(EL(i)));
    if (metal !== undefined) return { nH: null, charge: null, nHtotal: null, deprot: [], notes: [`contains ${EL(metal)}: metal complexes are not protonated by rules (load a MOL2 file)`], doubleBonds: [], metal: true };
    const nbOf = i => S.bonds[i].filter(j => S.atoms[j].el !== 'H');
    const bl = (i, j) => dist(pos(i), pos(j));
    const term = j => nbOf(j).length === 1;
    // planar unsaturated rings (5-7 atoms) inside the residue
    const ringSp2 = new Set();
    const adj = i => nbOf(i).filter(j => inR.has(j));
    const rings = []; const seen = new Set();
    for (const s0 of heavy) {
      const stack = [[s0, [s0]]];
      while (stack.length) {
        const [v, path] = stack.pop();
        if (path.length > 7) continue;
        for (const w of adj(v)) {
          if (w === s0 && path.length >= 5) { const key = [...path].sort((a, b) => a - b).join(','); if (!seen.has(key)) { seen.add(key); rings.push([...path]); } }
          else if (!path.includes(w) && w > s0) stack.push([w, path.concat(w)]);
        }
      }
    }
    for (const rg of rings) {
      const P0 = rg.map(pos); const c = P0.reduce((a, p) => add(a, p), [0, 0, 0]).map(x => x / P0.length);
      let nrm = [0, 0, 0]; for (let k = 0; k < P0.length; k++) nrm = add(nrm, cross(sub(P0[k], c), sub(P0[(k + 1) % P0.length], c)));
      if (norm(nrm) < 1e-6) continue; nrm = unit(nrm);
      const dev = Math.max(...P0.map(p => Math.abs(dot(sub(p, c), nrm))));
      let lsum = 0; for (let k = 0; k < rg.length; k++) lsum += bl(rg[k], rg[(k + 1) % rg.length]);
      if (dev < 0.12 && lsum / rg.length < 1.46) rg.forEach(i => ringSp2.add(i));
    }
    const cand = new Set(); const sp = new Set();
    // chain C with two heavy neighbours and trigonal angles (> 114 deg); used in pairs (C=C with poor deposited bond length)
    const trig2 = new Set(heavy.filter(i => EL(i) === 'C' && nbOf(i).length === 2 && !ringSp2.has(i) && angle(pos(nbOf(i)[0]), pos(i), pos(nbOf(i)[1])) > 114 && angle(pos(nbOf(i)[0]), pos(i), pos(nbOf(i)[1])) < 165));
    for (const i of heavy) {
      const e = EL(i), nb = nbOf(i), n = nb.length;
      const hyb = guessHyb(e, pos(i), nb.map(pos), nb.map(EL));
      if (hyb === 'sp' && (e === 'C' || e === 'N')) { sp.add(i); continue; }
      if (e === 'C') {
        const dmin = n ? Math.min(...nb.map(j => bl(i, j))) : 9;
        // chain C with two heavy neighbours: sp2 only with a double-bond length (a wide CH2 angle is not enough)
        const pair = trig2.has(i) && nb.some(j => trig2.has(j));
        if (ringSp2.has(i) || (n === 3 && hyb === 'sp2') || (n === 2 && hyb === 'sp2' && (dmin < 1.40 || pair)) || (n === 1 && hyb === 'sp2')) cand.add(i);
      }
      else if (e === 'N') { if (n === 2 && (ringSp2.has(i) || hyb === 'sp2')) cand.add(i); else if (n === 1 && EL(nb[0]) === 'C' && bl(i, nb[0]) < 1.32) cand.add(i); }
      else if ((e === 'O' || e === 'S') && n === 1 && EL(nb[0]) === 'C') { if (bl(i, nb[0]) < (e === 'O' ? 1.32 : 1.72)) cand.add(i); }
    }
    // Kekule structure by depth-first search over the sp2 carbons, one conjugated block at a time
    const partners = i => nbOf(i).filter(j => cand.has(j));
    const M = new Map(); let unmatchedC = 0, nodes = 0; const nodeCap = 200000;
    const comps = []; const vis = new Set();
    for (const s0 of cand) { if (vis.has(s0)) continue; const c = []; const st = [s0]; vis.add(s0); while (st.length) { const v = st.pop(); c.push(v); for (const w of partners(v)) if (!vis.has(w)) { vis.add(w); st.push(w); } } comps.push(c); }
    for (const comp of comps) {
      const mustC = comp.filter(i => EL(i) === 'C').sort((a, b) => partners(a).length - partners(b).length);
      let best = null; const mate = new Map(); let cn = 0;
      const score = () => { let o = 0, len = 0; for (const [i, j] of mate) if (i < j) { if (['O', 'S'].includes(EL(i)) || ['O', 'S'].includes(EL(j))) o++; len += bl(i, j); } return [mustC.filter(i => mate.has(i)).length, o, -len]; };
      const better = (a, b) => !b || a[0] > b[0] || (a[0] === b[0] && (a[1] > b[1] || (a[1] === b[1] && a[2] > b[2] + 1e-9)));
      const dfs = k => {
        if (++cn > nodeCap) return;
        if (k === mustC.length) {
          const extra = []; for (const i of comp) if (!mate.has(i) && EL(i) === 'N') for (const j of partners(i)) if (!mate.has(j) && !mate.has(i) && EL(j) !== 'C') { mate.set(i, j); mate.set(j, i); extra.push(i, j); }
          const sc = score(); if (better(sc, best && best.sc)) best = { sc, m: new Map(mate) };
          extra.forEach(i => mate.delete(i)); return;
        }
        const i = mustC[k];
        if (mate.has(i)) { dfs(k + 1); return; }
        for (const j of partners(i)) if (!mate.has(j)) { mate.set(i, j); mate.set(j, i); dfs(k + 1); mate.delete(i); mate.delete(j); }
        if (!best || best.sc[0] < mustC.length) dfs(k + 1); // leave i unmatched only while no full structure exists
      };
      dfs(0); nodes = Math.max(nodes, cn);
      if (best) { for (const [i, j] of best.m) M.set(i, j); unmatchedC += mustC.length - best.sc[0]; }
    }
    const nH = new Map(); const notes = [];
    for (const i of heavy) {
      const e = EL(i), n = nbOf(i).length, d = M.has(i) ? 1 : 0;
      let h = 0;
      if (e === 'C') h = sp.has(i) ? Math.max(0, 2 - n) : Math.max(0, 4 - n - d);
      else if (e === 'N') h = sp.has(i) ? 0 : (n >= 3 ? 0 : Math.max(0, 3 - n - d));
      else if (e === 'O' || e === 'S' || e === 'SE') {
        const nb = nbOf(i); const p = n === 1 ? nb[0] : -1;
        if (n >= 2 || (e === 'S' && n >= 2)) h = 0;
        else if (n === 1 && ['P', 'S', 'N', 'B', 'AS', 'SE'].includes(EL(p)) && e !== 'S') h = -1; // decided below
        else h = Math.max(0, 2 - n - d);
      }
      nH.set(i, h);
    }
    // terminal O on P / S / N
    const acidO = []; const deprot = [];
    for (const i of heavy) {
      const e = EL(i); if (!['P', 'S', 'N', 'AS', 'SE', 'B'].includes(e)) continue;
      const tO = nbOf(i).filter(j => EL(j) === 'O' && term(j)).sort((a, b) => bl(i, a) - bl(i, b));
      if (!tO.length) continue;
      if (e === 'N' || e === 'B') { tO.forEach(j => nH.set(j, e === 'B' ? 1 : 0)); continue; } // nitro / N-oxide: no H
      const nDouble = e === 'P' || e === 'AS' ? 1 : (nbOf(i).length >= 3 ? Math.min(2, tO.length) : 0);
      tO.forEach((j, k) => { if (k < nDouble) nH.set(j, 0); else { nH.set(j, 1); acidO.push({ c: i, o: j }); } });
    }
    for (const [i, h] of nH) if (h < 0) nH.set(i, 0);
    let charge = 0;
    if (mode !== 'zero') {
      for (const i of heavy) {
        const e = EL(i); const nb = nbOf(i);
        if (e === 'C' && nb.length >= 2 && nb.length <= 3) {
          const tO = nb.filter(j => EL(j) === 'O' && term(j));
          if (tO.length === 2) { const oh = tO.find(j => nH.get(j) === 1 && !M.has(j)); if (oh !== undefined) { nH.set(oh, 0); charge -= 1; deprot.push(oh); } }
        }
        if (e === 'P' || e === 'S') {
          const t = nb.filter(j => EL(j) === 'O' && term(j)).length;
          const want = e === 'P' ? (t >= 3 ? 2 : t === 2 ? 1 : 0) : (t >= 4 ? 2 : t === 3 ? 1 : 0);
          const ohs = acidO.filter(a => a.c === i).map(a => a.o).filter(j => nH.get(j) === 1);
          for (let k = 0; k < Math.min(want, ohs.length); k++) { nH.set(ohs[k], 0); charge -= 1; deprot.push(ohs[k]); }
        }
        if (e === 'N' && nb.length === 4 && !nb.some(j => EL(j) === 'O' && term(j))) charge += 1;
      }
    }
    if (unmatchedC) notes.push(`${unmatchedC} sp2 carbon(s) without a double-bond partner (unusual valence; check the ligand)`);
    if (nodes > nodeCap) notes.push('Kekule search stopped at the node limit');
    const total = [...nH.values()].reduce((a, b) => a + b, 0);
    return { nH, charge, nHtotal: total, deprot, notes, doubleBonds: [...M].filter(([i, j]) => i < j), sp2: new Set([...M.keys(), ...[...ringSp2].filter(i => cand.has(i))]) };
  }
  // FORMUL check: heavy-atom completeness and H count (expected H = FORMUL H + model charge - FORMUL charge)
  function formulCheck(S, R, nHmodel, charge) {
    const f = S.header && S.header.formul && S.header.formul[R.resName]; if (!f) return null;
    const cnt = {}; R.atoms.forEach(i => { const e = S.atoms[i].el; if (e !== 'H') cnt[e] = (cnt[e] || 0) + 1; });
    const heavyF = Object.entries(f).filter(([e]) => e !== 'H' && e !== '_charge').reduce((a, [, n]) => a + n, 0);
    const heavyM = Object.values(cnt).reduce((a, b) => a + b, 0);
    const complete = Object.entries(f).every(([e, n]) => e === 'H' || e === '_charge' || (cnt[e] || 0) === n);
    const inR = new Set(R.atoms); let links = 0; R.atoms.forEach(i => { if (S.atoms[i].el !== 'H') links += S.bonds[i].filter(j => !inR.has(j) && S.atoms[j].el !== 'H').length; });
    const expH = complete ? (f.H || 0) + charge - (f._charge || 0) - links : null;
    return { heavyF, heavyM, complete, links, expH, hOK: complete ? expH === nHmodel : null };
  }
  // Group rules for a ligand's net charge. Without hydrogens: carboxylate (two terminal O, C–O < 1.32 Å) −1,
  // phosphate −1/−2 (two/three terminal O), sulfonate −1, sulfate −2, quaternary N +1. With hydrogens present, an O carrying H
  // is not terminal and any N with four neighbours (incl. H) is ammonium (+1).
  function ligandChargeGuess(S, frame, R) {
    let q = 0; const hasH = R.atoms.some(i => S.atoms[i].el === 'H');
    for (const i of R.atoms) {
      const a = S.atoms[i]; const nb = S.bonds[i].filter(j => S.atoms[j].el !== 'H');
      const term = j => S.atoms[j].el === 'O' && S.bonds[j].length === 1;
      const termO = nb.filter(j => term(j) && dist(coordOf(frame, j), coordOf(frame, i)) < 1.32);
      if (a.el === 'C' && termO.length === 2) q -= 1;
      if (a.el === 'P') { const t = nb.filter(term); if (t.length >= 3) q -= 2; else if (t.length === 2) q -= 1; }
      if (a.el === 'S') { const t = nb.filter(term); if (t.length >= 4) q -= 2; else if (t.length === 3) q -= 1; }
      if (a.el === 'N' && (hasH ? S.bonds[i].length === 4 : nb.length === 4)) q += 1;
    }
    return q;
  }
  function linkPosition(pQ, pM, elQ, elM, P) {
    const xh = P.qm.xh[elQ] || 1.09; const dir = unit(sub(pM, pQ));
    if (P.qm.linkScheme === 'scaled') { const g = xh / ((COV[elQ] || 0.77) + (COV[elM] || 0.77)); return add(pQ, scl(sub(pM, pQ), g)); }
    return add(pQ, scl(dir, xh));
  }

  /* ------------------------------------------------------------------ */
  /*  T3a: QM model builder (cluster / QM-MM) with provenance manifest   */
  /* ------------------------------------------------------------------ */
  function buildQMModel(S, frameIdx, xbList, opts = {}) {
    const P = fullParams(opts.P); const frame = S.frames[frameIdx || 0]; const mode = P.qm.mode;
    const ffa = opts.ffa || assignForceField(S, P, { mol2: opts.mol2 });
    const FF = ffSet(P.ff.name);
    const donorResSet = new Set(xbList.map(h => h.donorRes)); const accResSet = new Set(xbList.map(h => h.accRes));
    const include = new Set(); const fragOf = new Map(); const partial = new Set();
    const resHasH = R => R.atoms.some(i => S.atoms[i].el === 'H');
    donorResSet.forEach(r => S.residues[r].atoms.forEach(i => { include.add(i); fragOf.set(i, 1); }));
    const BB = ['N', 'C', 'O', 'OXT', 'H', 'H1', 'H2', 'H3', 'HN', 'HT1', 'HT2', 'HT3', 'OT1', 'OT2'];
    accResSet.forEach(r => {
      if (donorResSet.has(r)) return;
      const R = S.residues[r]; const sc = mode === 'sidechain' && R.cls === 'protein' && (AA_ALIAS[R.resName] || R.resName) !== 'GLY';
      const needsBB = xbList.some(h => h.accRes === r && h.type === 'atom' && ['O', 'N', 'OXT'].includes(S.atoms[h.Y].name));
      R.atoms.forEach(i => { const nm = S.atoms[i].name; if (sc && !needsBB && BB.includes(nm)) { partial.add(r); return; } include.add(i); if (!fragOf.has(i)) fragOf.set(i, 2); });
    });
    if (P.qm.shell > 0) {
      const xs = xbList.map(h => coordOf(frame, h.X));
      for (const R of S.residues) {
        if (donorResSet.has(R.idx) || accResSet.has(R.idx)) continue;
        if (R.cls === 'ion' && !P.qm.shellIons) continue; if (R.cls === 'water' && !P.qm.shellWater) continue;
        if (R.atoms.some(i => xs.some(p => dist(p, coordOf(frame, i)) <= P.qm.shell))) R.atoms.forEach(i => { include.add(i); fragOf.set(i, 2); });
      }
    }
    const qm = []; const indexOf = new Map(); const qmRes = new Set([...include].map(i => S.atoms[i].res)); const warnings = ffa.warnings.filter(w => [...qmRes].some(r => w.includes(' ' + resId(S.residues[r]) + ':')));
    const droppedH = new Set(); // file H replaced by rebuilt H: still QM atoms for the embedding (never MM charges)
    if (P.qm.forceH) [...include].forEach(i => { if (S.atoms[i].el === 'H') { include.delete(i); fragOf.delete(i); droppedH.add(i); } }); // re-protonate: drop existing H
    // ligands: one heavy-atom structure gives both the hydrogens and the charge (ligandProtonation);
    // hydrogens present in the file are kept only when their number matches that structure
    const ligInfo = {}; const ligRebuild = new Set();
    [...new Set([...include].map(i => S.atoms[i].res))].forEach(r => {
      const R = S.residues[r]; if (R.cls !== 'ligand') return;
      const lp = ligandProtonation(S, frame, R, P.qm.ligandCharge === 'zero' ? 'zero' : 'rules');
      const fileH = R.atoms.filter(i => S.atoms[i].el === 'H').length; const id = `${R.resName} ${resId(R)}`;
      const L = { res: resId(R), resName: R.resName, fileH, rebuilt: false, notes: lp.notes.slice() };
      if (lp.metal) { L.notes.forEach(n => warnings.push(`${id}: ${n}; old per-atom rules used.`)); ligInfo[r] = { lp: null, L }; return; }
      if (fileH && !P.qm.forceH && fileH !== lp.nHtotal) {
        R.atoms.forEach(i => { if (S.atoms[i].el === 'H') { include.delete(i); fragOf.delete(i); droppedH.add(i); } }); ligRebuild.add(r); L.rebuilt = true;
        warnings.push(`${id}: the file has ${fileH} H, but its heavy-atom structure needs ${lp.nHtotal} H at charge ${lp.charge}; its hydrogens were rebuilt.`);
      }
      L.modelH = (fileH && !P.qm.forceH && !L.rebuilt) ? fileH : lp.nHtotal; L.charge = lp.charge;
      const fc = formulCheck(S, R, L.modelH, lp.charge); L.formul = fc;
      if (fc && !fc.complete) warnings.push(`${id}: ${fc.heavyM} of ${fc.heavyF} heavy atoms of its PDB formula (FORMUL) are present; the missing atoms are not modelled.`);
      else if (fc && !fc.hOK) warnings.push(`${id}: ${L.modelH} H at charge ${lp.charge}, but the PDB formula implies ${fc.expH} H; check the ligand protonation or load a MOL2 file.`);
      lp.notes.forEach(n => warnings.push(`${id}: ${n}.`));
      ligInfo[r] = { lp, L };
    });
    const order = [...include].sort((a, b) => (fragOf.get(a) - fragOf.get(b)) || a - b);
    for (const i of order) { indexOf.set(i, qm.length); const a = S.atoms[i]; const p = coordOf(frame, i);
      qm.push({ el: a.el, x: p[0], y: p[1], z: p[2], frag: fragOf.get(i), src: { i, chain: a.chain, resName: a.resName, resSeq: a.resSeq, iCode: a.iCode, name: a.name, serial: a.serial } }); }
    let addedH = 0, linkH = 0, completedH = 0; const keptH = [];
    for (const i of order) {
      const a = S.atoms[i]; if (a.el === 'H') continue; const R = S.residues[a.res]; const info = ffa.residues[R.idx];
      const pc = coordOf(frame, i); const nbAll = S.bonds[i].filter(j => S.atoms[j].el !== 'H');
      const protonate = P.qm.forceH || ligRebuild.has(R.idx) || !resHasH(R); const LP = ligInfo[R.idx] && ligInfo[R.idx].lp;
      if (!protonate && !keptH.includes(R.idx)) keptH.push(R.idx);
      // Residues that already carry H: N/O/S hydrogens are taken as given; carbons that lack their H
      // (polar-hydrogen-only files) are completed so that the QM model is never under-protonated.
      const existingH = S.bonds[i].filter(j => S.atoms[j].el === 'H');
      if (protonate || a.el === 'C') {
        let nH = 0, hyb = 'sp3'; const tpl = info && info.tpl && FF.residues[info.tpl];
        if (R.cls === 'water') { nH = 2; hyb = 'sp3'; }
        else if (R.cls === 'ion') nH = 0;
        else if (tpl && (R.cls === 'protein' || R.cls === 'nucleic')) {
          let n = a.name; if (R.resName === 'MSE' && n === 'SE') n = 'SD'; if (!tpl.atoms[n] && HEAVY_ALIAS[n]) n = HEAVY_ALIAS[n]; if (!tpl.atoms[n] && R.resName === 'ILE' && n === 'CD') n = 'CD1';
          const t = tpl.atoms[n]; nH = t ? t[2].length : 0;
          const rn = AA_ALIAS[R.resName] || R.resName; const sp2 = SP2_ATOMS._BB.includes(a.name) || (SP2_ATOMS[rn] || []).includes(a.name);
          hyb = R.cls === 'nucleic' ? guessHyb(a.el, pc, nbAll.map(j => coordOf(frame, j)), nbAll.map(j => S.atoms[j].el)) : (sp2 ? 'sp2' : 'sp3');
          if (a.name === 'N' && nH >= 2) hyb = 'sp3';
          if (rn === 'TYR' && a.name === 'OH') hyb = 'sp2';
        } else if (LP) {
          const nbPos = nbAll.map(j => coordOf(frame, j)); const nbEl = nbAll.map(j => S.atoms[j].el);
          nH = LP.nH.get(i) || 0; hyb = LP.sp2.has(i) ? 'sp2' : guessHyb(a.el, pc, nbPos, nbEl);
          if (a.el === 'C' && nH + nbAll.length >= 4) hyb = 'sp3'; // saturated carbon: tetrahedral H
        } else {
          const nbPos = nbAll.map(j => coordOf(frame, j)); const nbEl = nbAll.map(j => S.atoms[j].el);
          hyb = guessHyb(a.el, pc, nbPos, nbEl); nH = heuristicHCount(a.el, hyb, nbAll.length, pc, nbPos);
        }
        if (!protonate) nH -= existingH.filter(j => include.has(j)).length;
        if (nH > 0) {
          const nbPos = nbAll.map(j => coordOf(frame, j)).concat(protonate ? [] : existingH.filter(j => include.has(j)).map(j => coordOf(frame, j)));
          let ref = null; if (nbPos.length === 1 && nbAll.length === 1) { const k = S.bonds[nbAll[0]].find(j => j !== i && S.atoms[j].el !== 'H'); if (k !== undefined) ref = coordOf(frame, k); }
          placeH(pc, nbPos, nH, hyb, ref, P.qm.xh[a.el] || 1.0).forEach(h => { qm.push({ el: 'H', x: h[0], y: h[1], z: h[2], frag: fragOf.get(i), added: true, parent: indexOf.get(i) }); addedH++; if (!protonate) completedH++; });
        }
      }
      for (const j of nbAll) if (!include.has(j)) {
        const h = linkPosition(pc, coordOf(frame, j), a.el, S.atoms[j].el, P);
        qm.push({ el: 'H', x: h[0], y: h[1], z: h[2], frag: fragOf.get(i), link: true, parent: indexOf.get(i), replaces: { i: j, label: atomLabel(S, j) } }); linkH++;
      }
    }
    if (keptH.length) warnings.push(`${keptH.length} residue(s) already contain hydrogens (${keptH.slice(0, 5).map(r => S.residues[r].resName + ' ' + resId(S.residues[r])).join(', ')}${keptH.length > 5 ? '…' : ''}): their N/O/S hydrogens are kept as given${completedH ? `; ${completedH} missing C–H hydrogen(s) were added` : ' and no H is added'}.`);
    // ---- fragment charges from force-field formal charges
    const fragCharge = { 1: 0, 2: 0 }; const chargeLog = []; const qmResCharge = {};
    const seenRes = [...new Set(order.map(i => S.atoms[i].res))];
    seenRes.forEach(r => {
      const R = S.residues[r]; const info = ffa.residues[r]; const fr = fragOf.get(R.atoms.find(i => include.has(i)));
      let c = 0, how = '';
      if (R.cls === 'ligand') {
        const m = opts.mol2 && opts.mol2[R.resName];
        if (P.qm.ligandCharge === 'zero') { c = 0; how = 'set to zero'; }
        else if (P.qm.ligandCharge === 'auto' && m) { c = Math.round(m._total); how = 'MOL2 total charge'; }
        else if (info && info.source === 'PQR' && !(ligInfo[r] && ligInfo[r].L.rebuilt)) { c = info.formal; how = 'PQR charges'; }
        else if (ligInfo[r] && ligInfo[r].lp) { c = ligInfo[r].lp.charge; how = 'heavy-atom structure + group rules (COO⁻, phosphate, sulfonate/sulfate, R₄N⁺)'; }
        else { c = ligandChargeGuess(S, frame, R); how = 'group rules (carboxylate, phosphate, sulfonate/sulfate, R₄N⁺)'; }
      } else if (info && info.formal !== null && info.formal !== undefined) {
        c = info.formal; how = `${info.tpl || info.resName} template`;
        if (partial.has(r)) { // side chain only: remove terminal backbone charges
          let t = 0; if (info.tpl && /^N/.test(info.tpl) && info.tpl.length === 4) t += 1; if (info.tpl && /^C/.test(info.tpl) && info.tpl.length === 4 && !['CYS', 'CYX', 'CYM'].includes(info.tpl)) t -= 1;
          c -= t; how += ' (side chain only)';
        }
      } else { c = 0; how = 'unknown (0)'; warnings.push(`${R.resName} ${resId(R)} in the QM region has no charge information; 0 assumed.`); }
      fragCharge[fr] = (fragCharge[fr] || 0) + c; qmResCharge[r] = c; chargeLog.push({ res: resId(R), resName: R.resName, frag: fr, tpl: info ? info.tpl : null, charge: c, how });
    });
    if (opts.chargeOverride) Object.assign(fragCharge, opts.chargeOverride);
    const totalQ = (fragCharge[1] || 0) + (fragCharge[2] || 0);
    const electrons = qm.reduce((s, a) => s + (Z[a.el] || 0), 0) - totalQ;
    const fragElectrons = f => qm.filter(a => a.frag === f).reduce((s, a) => s + (Z[a.el] || 0), 0) - (fragCharge[f] || 0);
    if (electrons % 2 !== 0) warnings.push(`Odd electron count (${electrons}) for total charge ${totalQ}: check protonation/charges or set multiplicity 2.`);
    [1, 2].forEach(f => { const e = fragElectrons(f); if (e % 2 !== 0) warnings.push(`Fragment ${f} has an odd electron count (${e}); adjust its charge before ALMO-EDA.`); });
    const xbs = xbList.map((h, k) => ({ id: h.cid || ('XB' + (k + 1)), type: h.type, X: indexOf.get(h.X) + 1, A: indexOf.get(h.A) + 1,
      Y: h.type === 'atom' ? indexOf.get(h.Y) + 1 : null, ring: h.type === 'pi' ? h.ring.map(i => indexOf.get(i) + 1) : null,
      label: `${h.donorLabel}@${h.Xname} ··· ${h.accLabel}@${h.accName}`, d: h.d, theta1: h.theta1, halogen: h.halogen }));
    const manifest = { tool: 'BondX', version: VERSION, created: new Date().toISOString(), source: S.header.idCode || opts.sourceName || '', frame: frameIdx || 0, mode,
      shell: P.qm.shell, nAtoms: qm.length, addedH, linkH, forceField: P.ff.name, fragCharge, totalCharge: totalQ, multiplicity: opts.multiplicity || (electrons % 2 ? 2 : 1),
      residueCharges: chargeLog, ligands: Object.values(ligInfo).map(x => x.L), xbs,
      atoms: qm.map((a, k) => ({ n: k + 1, el: a.el, frag: a.frag, x: +a.x.toFixed(5), y: +a.y.toFixed(5), z: +a.z.toFixed(5),
        src: a.src ? `${a.src.resName}${a.src.resSeq}${a.src.iCode || ''}:${a.src.chain}@${a.src.name}` : (a.link ? `link-H(${a.parent + 1})→${a.replaces.label}` : `H(${a.parent + 1})`) })),
      warnings, parameters: { qm: P.qm, ff: Object.assign({}, P.ff, { overrides: P.ff.overrides }) }, assumptions: assumptions('T3', Object.assign({}, P, { qm: Object.assign({}, P.qm, { env: 'cluster' }) })) };
    manifest.id = 'BX-' + hashString(JSON.stringify(manifest.atoms) + JSON.stringify(xbs)).slice(0, 10);
    // include = structure atoms of the QM region (used by the embedding); replaced file H stay in it
    return { qm, manifest, include: new Set([...include, ...droppedH]), indexOf, fragOf, qmResCharge, ffa };
  }
  function hashString(s) { let h1 = 0x811c9dc5, h2 = 0x01000193; for (let i = 0; i < s.length; i++) { const c = s.charCodeAt(i); h1 = Math.imul(h1 ^ c, 16777619) >>> 0; h2 = Math.imul(h2 ^ c, 2246822519) >>> 0; } return h1.toString(16).padStart(8, '0') + h2.toString(16).padStart(8, '0'); }

  /* ---------------- Embedding point charges (QM/MM) ------------------ */
  function embeddingCharges(S, frameIdx, include, opts = {}) {
    const P = fullParams(opts.P); const frame = S.frames[frameIdx || 0]; const R = P.mm.radius; const ffa = opts.ffa || assignForceField(S, P, { mol2: opts.mol2 });
    const qmIdx = [...include]; const grid = makeGrid(frame, qmIdx, Math.max(4, R / 2));
    const near = i => { const p = coordOf(frame, i); let hit = false; gridQuery(grid, p, R, j => { if (!hit && dist(p, coordOf(frame, j)) <= R) hit = true; }); return hit; };
    // MM atoms
    const mmSet = new Set();
    if (P.mm.selection === 'residue') S.residues.forEach(Rs => { if (Rs.atoms.some(i => !include.has(i) && near(i))) Rs.atoms.forEach(i => { if (!include.has(i)) mmSet.add(i); }); });
    else S.atoms.forEach(a => { if (!include.has(a.i) && near(a.i)) mmSet.add(a.i); });
    const q = new Map(); mmSet.forEach(i => q.set(i, ffa.q[i]));
    const notes = [];
    // cut residues: residues with atoms on both sides -> MM part gets an integer charge (formal - QM-part formal)
    const cut = []; const qmFormal = opts.qmResidueCharges || {}; // res idx -> formal charge carried by QM part
    S.residues.forEach(Rs => { const inQ = Rs.atoms.filter(i => include.has(i)), inM = Rs.atoms.filter(i => mmSet.has(i)); if (inQ.length && inM.length) cut.push({ Rs, inQ, inM }); });
    if (P.mm.neutralizeCut) cut.forEach(({ Rs, inM }) => {
      const info = ffa.residues[Rs.idx]; if (!info || info.formal === null) return;
      const qmPart = qmFormal[Rs.idx] !== undefined ? qmFormal[Rs.idx] : Math.round(Rs.atoms.filter(i => include.has(i)).reduce((s, i) => s + ffa.q[i], 0));
      const target = info.formal - qmPart; const cur = inM.reduce((s, i) => s + q.get(i), 0); const dq = (target - cur) / inM.length;
      inM.forEach(i => q.set(i, q.get(i) + dq)); notes.push(`${resId(Rs)}: MM part of the cut residue shifted by ${(target - cur).toFixed(4)} e to carry ${target}.`);
    });
    // boundary scheme
    const M1 = new Set(); include.forEach(i => S.bonds[i].forEach(j => { if (mmSet.has(j)) M1.add(j); }));
    const nextShell = (prev, excl) => { const s = new Set(); prev.forEach(i => S.bonds[i].forEach(j => { if (mmSet.has(j) && !excl.has(j)) s.add(j); })); return s; };
    const M2 = nextShell(M1, new Set([...M1])); const M3 = nextShell(M2, new Set([...M1, ...M2]));
    const extraPts = []; const sch = P.mm.boundary;
    if (sch === 'Z1' || sch === 'Z2' || sch === 'Z3') { const del = [...M1, ...(sch !== 'Z1' ? M2 : []), ...(sch === 'Z3' ? M3 : [])].reduce((t, i) => t + q.get(i), 0);
      if (Math.abs(del) > 1e-6) notes.push(`${sch} boundary deleted ${del.toFixed(4)} e of MM charge near the QM region, so the net MM charge is no longer an integer (inherent to ${sch}; RC conserves it).`);
      M1.forEach(i => q.set(i, 0)); if (sch !== 'Z1') M2.forEach(i => q.set(i, 0)); if (sch === 'Z3') M3.forEach(i => q.set(i, 0)); }
    else if (sch === 'RC') {
      M1.forEach(i => { const qi = q.get(i); const m2 = S.bonds[i].filter(j => mmSet.has(j) && !M1.has(j)); q.set(i, 0); if (!m2.length) { if (qi) notes.push(`${atomLabel(S, i)}: M1 atom without M2 neighbours, its charge (${qi.toFixed(3)}) was deleted.`); return; }
        const pi = coordOf(frame, i); m2.forEach(j => { const pm = scl(add(pi, coordOf(frame, j)), 0.5); extraPts.push({ x: pm[0], y: pm[1], z: pm[2], q: qi / m2.length, label: `RC ${atomLabel(S, i)}–${S.atoms[j].name}` }); }); });
    }
    const out = []; mmSet.forEach(i => { const v = q.get(i); if (Math.abs(v) > 1e-8) { const p = coordOf(frame, i); const a = S.atoms[i]; out.push({ x: p[0], y: p[1], z: p[2], q: v, i, label: atomLabel(S, i), type: ffa.type[i] || '', src: ffa.residues[a.res].source }); } });
    extraPts.forEach(e => out.push(e));
    const unknown = [...new Set([...mmSet].map(i => S.atoms[i].res))].map(r => ffa.residues[r]).filter(x => x && x.source === 'none');
    if (unknown.length) notes.push(`${unknown.length} MM residue(s) without charges (set to 0): ${[...new Set(unknown.map(u => u.resName))].join(', ')}. Load a MOL2/PQR file to include them.`);
    const formal = [...new Set([...mmSet].map(i => S.atoms[i].res))].map(r => ffa.residues[r]).filter(x => x && x.source === 'formal');
    if (formal.length) notes.push(`${formal.length} MM residue(s) without template or MOL2 carry their formal charge spread over heavy atoms: ${[...new Set(formal.map(u => u.resName + ' ' + (u.formal > 0 ? '+' : '') + u.formal))].join(', ')}. Load a MOL2 file for proper partial charges.`);
    const missingH = P.ff.foldH && !ffa.usedPQR && [...mmSet].some(i => (ffa.residues[S.atoms[i].res] || {}).foldedH > 0);
    return { charges: out, source: ffa.usedPQR ? 'PQR' : ffa.ff, water: ffa.water, net: out.reduce((s, c) => s + c.q, 0), nMM: mmSet.size, M1: M1.size, M2: M2.size, cut: cut.length,
      boundary: sch, radius: R, selection: P.mm.selection, notes, foldedH: missingH };
  }

  /* ---------------- QM input writers ---------------------------------- */
  const asc = s => String(s).replace(/···/g, '...').replace(/[—–]/g, '-').replace(/→/g, '->').replace(/π/g, 'pi').replace(/[^\x20-\x7E]/g, '?');
  const fmtXYZ = a => `${a.el.charAt(0) + a.el.slice(1).toLowerCase()}`.padEnd(3) + a.x.toFixed(6).padStart(13) + a.y.toFixed(6).padStart(13) + a.z.toFixed(6).padStart(13);
  function heavyECP(qm) { return qm.some(a => (Z[a.el] || 0) > 36); }
  // Dispersion correction to add for a method: built-in (wB97X-D, -V, D3/D4 names), none (HF, MP2, CC, Minnesota),
  // D3(BJ) for other hybrids. o.disp ('none' | 'bj' | 'zero') overrides the automatic choice.
  function dispKind(method, explicit) {
    if (explicit) return explicit;
    if (/XD|-D3|-D4|D3BJ|-V$|^HF$|MP2|CCSD|DLPNO|^M0[5-8]|^M11|^MN1|^SOGGA|^B97-?D$/i.test(method)) return 'none';
    return 'bj';
  }
  function writeGaussian(model, o, pcs) {
    const m = model.manifest; const q = m.totalCharge, mult = m.multiplicity;
    const dk = dispKind(o.method, o.disp); const disp = dk === 'bj' ? ' empiricaldispersion=gd3bj' : dk === 'zero' ? ' empiricaldispersion=gd3' : '';
    const solv = (!pcs && o.solvent) ? ' scrf=(cpcm,solvent=water,read)' : '';
    const L = [`%nprocshared=${o.nproc}`, `%mem=${o.mem}`, `%chk=${o.job}.chk`,
      `#p ${o.method}/${o.basis}${disp} pop=(nbo${o.nbo7 ? '7' : ''}) nosymm scf=(tight,xqc) int=ultrafine${pcs ? ' charge' : ''}${solv}`, '',
      asc(`BondX ${m.id} | NBO XB gate | ${m.xbs.map(x => x.label).join(' ; ')}`).slice(0, 78), '', `${q} ${mult}`];
    model.qm.forEach(a => L.push(fmtXYZ(a))); L.push('');
    if (pcs) { pcs.charges.forEach(c => L.push(`${c.x.toFixed(6).padStart(13)}${c.y.toFixed(6).padStart(13)}${c.z.toFixed(6).padStart(13)}${c.q.toFixed(6).padStart(12)}`)); L.push(''); }
    if (solv) L.push(`Eps=${o.eps}`, '');
    L.push('');
    return L.join('\n');
  }
  function writeORCA(model, o, pcs) {
    const m = model.manifest; const L = [`# BondX ${m.id} - NBO XB gate`, asc(`# ${m.xbs.map(x => x.label).join(' ; ')}`),
      '# CAUTION: ORCA has no built-in NBO. The NBO keyword calls NBO 7, a licensed program.',
      '# Set the NBOEXE environment variable to the NBO 7 executable before running ORCA.',
      '# Without NBO 7 the SCF runs but no NBO analysis is printed.',
      '# The ORCA + NBO 7 route is not tested in this version of BondX.'];
    const dk = dispKind(o.method, o.disp); const disp = dk === 'bj' ? ' D3BJ' : dk === 'zero' ? ' D3ZERO' : '';
    L.push(`! ${o.method} ${o.basis}${disp} RIJCOSX def2/J TightSCF NBO${(!pcs && o.solvent) ? ' CPCM(Water)' : ''}`);
    L.push(`%pal nprocs ${o.nproc} end`, `%maxcore ${o.maxcore}`);
    if (!pcs && o.solvent) L.push('%cpcm', `  epsilon ${o.eps}`, 'end');
    if (pcs) L.push(`%pointcharges "${o.job}.pc"`);
    L.push('', `* xyz ${m.totalCharge} ${m.multiplicity}`); model.qm.forEach(a => L.push(fmtXYZ(a))); L.push('*', '');
    const pc = pcs ? [String(pcs.charges.length), ...pcs.charges.map(c => `${c.q.toFixed(6)} ${c.x.toFixed(6)} ${c.y.toFixed(6)} ${c.z.toFixed(6)}`)].join('\n') + '\n' : null;
    return { inp: L.join('\n'), pc };
  }
  function qchemRem(o, extra, qm) {
    const r = [`   METHOD          ${o.method}`, `   BASIS           ${o.basis}`];
    if (o.disp === 'bj') r.push('   DFT_D           D3_BJ'); else if (o.disp === 'zero') r.push('   DFT_D           D3_ZERO');
    if (heavyECP(qm) && /def2/i.test(o.basis)) r.push('   ECP             def2-ECP');
    r.push('   SCF_CONVERGENCE 8', '   THRESH          14', '   SYM_IGNORE      TRUE', '   XC_GRID         000099000590', '   MEM_TOTAL       ' + (o.memMB || 16000));
    if (/-V$/i.test(o.method)) r.push('   NL_GRID         1');
    return r.concat(extra);
  }
  function writeQChemNBO(model, o, pcs) {
    const m = model.manifest; const L = [`! BondX ${m.id} - NBO XB gate`, asc(`! ${m.xbs.map(x => x.label).join(' ; ')}`)];
    if (QCHEM_G_BASIS.test(o.basis)) L.push('! ' + QCHEM_NBO_G_CAUTION);
    L.push('$molecule', `${m.totalCharge} ${m.multiplicity}`);
    model.qm.forEach(a => L.push(fmtXYZ(a))); L.push('$end', '');
    const rem = ['$rem', '   JOBTYPE         sp', ...qchemRem(o, ['   NBO             1'], model.qm)];
    if (!pcs && o.solvent) rem.push('   SOLVENT_METHOD  PCM');
    rem.push('$end'); L.push(...rem);
    // NBO 1 alone prints the NBOs but not the second-order table in Q-Chem 6.2 (NBO 5.0): request it explicitly
    L.push('', '$nbo', '   E2PERT', '$end');
    if (pcs) { L.push('', '$external_charges'); pcs.charges.forEach(c => L.push(`${c.x.toFixed(6)} ${c.y.toFixed(6)} ${c.z.toFixed(6)} ${c.q.toFixed(6)}`)); L.push('$end'); }
    else if (o.solvent) L.push('', '$solvent', `   Dielectric ${o.eps}`, '$end');
    return L.join('\n') + '\n';
  }
  function writeQChemEDA(model, o) {
    const m = model.manifest; const fq = m.fragCharge;
    const e1 = model.qm.filter(a => a.frag === 1).reduce((s, a) => s + (Z[a.el] || 0), 0) - (fq[1] || 0);
    const e2 = model.qm.filter(a => a.frag !== 1).reduce((s, a) => s + (Z[a.el] || 0), 0) - (fq[2] || 0);
    const L = [`! BondX ${m.id} - ALMO-EDA (EDA2) on NBO-verified XB(s)`, asc(`! ${m.xbs.map(x => x.label).join(' ; ')}`),
      asc(`! Fragment 1 = halogen-bond donor (${[...new Set(m.atoms.filter(a => a.frag === 1 && !a.src.startsWith('H') && !a.src.startsWith('link')).map(a => a.src.split('@')[0]))].join(', ')})`),
      `! Fragment 2 = halogen-bond acceptor environment`, '$molecule', `${m.totalCharge} ${m.multiplicity}`, '--', `${fq[1] || 0} ${e1 % 2 ? 2 : 1}`];
    model.qm.filter(a => a.frag === 1).forEach(a => L.push(fmtXYZ(a)));
    L.push('--', `${fq[2] || 0} ${e2 % 2 ? 2 : 1}`); model.qm.filter(a => a.frag !== 1).forEach(a => L.push(fmtXYZ(a))); L.push('$end', '');
    const extra = ['   EDA2            ' + (o.eda2 || 1), '   EDA_BSSE        ' + (o.bsse ? 'TRUE' : 'FALSE')];
    if (o.covp) extra.push('   EDA_COVP        TRUE'); if (o.solvent) extra.push('   SOLVENT_METHOD  PCM'); if (e1 % 2 || e2 % 2) extra.push('   UNRESTRICTED    TRUE');
    L.push('$rem', '   JOBTYPE         eda', ...qchemRem(o, extra, model.qm), '$end');
    if (o.solvent) L.push('', '$solvent', `   SolventDielectric ${o.eps}`, '$end');
    return L.join('\n') + '\n';
  }
  function writeXYZ(model) { return `${model.qm.length}\nBondX ${model.manifest.id}\n` + model.qm.map(fmtXYZ).join('\n') + '\n'; }
  function chargesCSV(pcs) { return 'x,y,z,q,atom,type,source\n' + pcs.charges.map(c => [c.x.toFixed(4), c.y.toFixed(4), c.z.toFixed(4), c.q.toFixed(6), JSON.stringify(c.label || ''), c.type || '', c.src || 'RC'].join(',')).join('\n') + '\n'; }

  /* ------------------------------------------------------------------ */
  /*  T3b: NBO output parsing and the XB gate                            */
  /* ------------------------------------------------------------------ */
  const NBO_SIDE = String.raw`(\d+|\*+)\.\s*(BD\*?|LP\*?|CR|RY\*?|LV|3C\*?)\s*\(\s*(\d+)\s*\)\s*([A-Z][a-z]?)\s*(\d+)(?:\s*-\s*([A-Z][a-z]?)\s*(\d+))?(?:\s*-\s*([A-Z][a-z]?)\s*(\d+))?`;
  const NBO_LINE = new RegExp('^\\s*' + NBO_SIDE + '\\s*\\/\\s*' + NBO_SIDE + '\\s+(\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)');
  function parseNBOSide(m, o) { const atoms = []; for (let k = 0; k < 3; k++) { if (m[o + 3 + 2 * k]) atoms.push({ el: m[o + 3 + 2 * k].toUpperCase(), n: parseInt(m[o + 4 + 2 * k]) }); } return { nbo: /^\*+$/.test(m[o]) ? null : parseInt(m[o]), type: m[o + 1].replace(/\s/g, ''), idx: parseInt(m[o + 2]), atoms }; }
  function parseNBO(text) {
    const lines = text.split(/\r?\n/);
    const program = /Entering Gaussian System|Gaussian, Inc/.test(text) ? 'Gaussian' : /O\s+R\s+C\s+A/.test(text) || /ORCA TERMINATED|\* O   R   C   A \*/.test(text) ? 'ORCA' : /Q-Chem/.test(text) ? 'Q-Chem' : 'NBO';
    const terminated = /Normal termination of Gaussian|ORCA TERMINATED NORMALLY|Thank you very much for using Q-Chem|Total job time/.test(text);
    const secStarts = []; lines.forEach((l, i) => { if (/SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS/i.test(l)) secStarts.push(i); });
    if (!secStarts.length) return { program, terminated, entries: [], error: 'No "SECOND ORDER PERTURBATION THEORY ANALYSIS" section found. Ensure NBO was requested (Gaussian pop=nbo, ORCA ! NBO with NBO 7, Q-Chem NBO 1 with $nbo E2PERT).' };
    const spinOf = i => { for (let k = i; k >= Math.max(0, i - 40); k--) { if (/Alpha spin orbitals/i.test(lines[k])) return 'alpha'; if (/Beta\s+spin orbitals/i.test(lines[k])) return 'beta'; } return 'closed'; };
    let use = [secStarts[secStarts.length - 1]];
    if (spinOf(use[0]) === 'beta' && secStarts.length >= 2 && spinOf(secStarts[secStarts.length - 2]) === 'alpha') use = secStarts.slice(-2);
    const entries = []; let threshold = null;
    for (const s of use) {
      const spin = spinOf(s);
      for (let i = s + 1; i < lines.length; i++) {
        const l = lines[i];
        if (/Threshold for printing:\s*([\d.]+)/.test(l)) threshold = parseFloat(RegExp.$1);
        if (/NATURAL BOND ORBITALS \(Summary\)|NBO\s+.*Summary/i.test(l) && i > s + 5) break;
        const m = NBO_LINE.exec(l);
        if (m) entries.push({ donor: parseNBOSide(m, 1), acceptor: parseNBOSide(m, 10), e2: parseFloat(m[19]), dE: parseFloat(m[20]), F: parseFloat(m[21]), spin });
      }
    }
    const npa = {}; let iN = -1; lines.forEach((l, i) => { if (/Summary of Natural Population Analysis/i.test(l)) iN = i; });
    if (iN >= 0) for (let i = iN + 1; i < Math.min(lines.length, iN + 4000); i++) {
      const m = /^\s*([A-Z][a-z]?)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/.exec(lines[i]);
      if (m) npa[parseInt(m[2])] = parseFloat(m[3]); else if (/^\s*=+\s*$/.test(lines[i]) && Object.keys(npa).length) break; else if (/\* Total \*/.test(lines[i])) break;
    }
    return { program, terminated, entries, threshold, npa, geometry: parseOutputGeometry(text, program) };
  }
  function parseOutputGeometry(text, program) {
    const lines = text.split(/\r?\n/); let geo = null;
    if (program === 'Gaussian') { let st = -1; lines.forEach((l, i) => { if (/Input orientation:|Standard orientation:/.test(l)) st = i; });
      if (st >= 0) { geo = []; for (let j = st + 5; j < lines.length; j++) { if (/^\s*-{5,}/.test(lines[j])) break; const t = lines[j].trim().split(/\s+/); geo.push({ el: ZSYM[+t[1]] || 'X', x: +t[3], y: +t[4], z: +t[5] }); } } }
    else if (program === 'ORCA') { let st = -1; lines.forEach((l, i) => { if (/CARTESIAN COORDINATES \(ANGSTROEM\)/.test(l)) st = i; });
      if (st >= 0) { geo = []; for (let j = st + 2; j < lines.length; j++) { const t = lines[j].trim().split(/\s+/); if (t.length < 4 || isNaN(+t[1])) break; geo.push({ el: t[0].toUpperCase(), x: +t[1], y: +t[2], z: +t[3] }); } } }
    else if (program === 'Q-Chem') { let st = -1; lines.forEach((l, i) => { if (/Standard Nuclear Orientation/.test(l)) st = i; });
      if (st >= 0) { geo = []; for (let j = st + 3; j < lines.length; j++) { if (/^\s*-{5,}/.test(lines[j])) break; const t = lines[j].trim().split(/\s+/); geo.push({ el: t[1].toUpperCase(), x: +t[2], y: +t[3], z: +t[4] }); } } }
    return geo;
  }
  function verifyGeometry(manifest, geo, tol = DEFAULTS.nbo.geomTol) {
    if (!geo || !geo.length) return { checked: false, note: 'No geometry block found in output; atom-order verification skipped.' };
    if (geo.length !== manifest.atoms.length) return { checked: true, ok: false, elMis: Math.abs(geo.length - manifest.atoms.length), note: `Atom count mismatch: output ${geo.length} vs manifest ${manifest.atoms.length}.` };
    let maxDev = 0, elMis = 0;
    manifest.atoms.forEach((a, k) => { const g = geo[k]; if (g.el !== a.el) elMis++; maxDev = Math.max(maxDev, dist([a.x, a.y, a.z], [g.x, g.y, g.z])); });
    return { checked: true, ok: elMis === 0 && maxDev < tol, maxDev, elMis, note: elMis ? `${elMis} element mismatches — output does not belong to this manifest.` : (maxDev < tol ? `Atom order verified (max deviation ${maxDev.toFixed(4)} Å).` : `Coordinates differ by up to ${maxDev.toFixed(3)} Å (> ${tol} Å; re-oriented or optimised?). Element order matches.`) };
  }
  function nboGate(manifest, nbo, opts = {}) {
    const thr = opts.threshold ?? DEFAULTS.nbo.threshold; const useRY = opts.includeRY ?? DEFAULTS.nbo.includeRY;
    return manifest.xbs.map(xb => {
      const isA = s => s.atoms.length === 2 && s.atoms.some(a => a.n === xb.X) && s.atoms.some(a => a.n === xb.A);
      const ringSet = new Set(xb.ring || []);
      const contrib = { lp_sigma: 0, lp_ry: 0, pi_sigma: 0, other: 0 }; const rows = [];
      // Ionic Lewis structure: NBO may describe A-X as A(-) + X(+) with a vacant lone pair LP*(X) instead of sigma*(A-X).
      // Then the donation into LP*(X) takes the role of the n -> sigma* term.
      const hasSigmaStar = nbo.entries.some(e => e.acceptor.type === 'BD*' && isA(e.acceptor));
      const isLPstarX = s => !hasSigmaStar && s.type === 'LP*' && s.atoms.length === 1 && s.atoms[0].n === xb.X;
      for (const e of nbo.entries) {
        const d = e.donor, a = e.acceptor;
        const toSigma = (a.type === 'BD*' && isA(a)) || isLPstarX(a); const toRyX = a.type.startsWith('RY') && a.atoms.length === 1 && a.atoms[0].n === xb.X;
        if (!toSigma && !toRyX) continue;
        let tag = null;
        if (xb.type === 'atom' && d.type === 'LP' && d.atoms[0].n === xb.Y) tag = toSigma ? 'lp_sigma' : 'lp_ry';
        else if (xb.type === 'pi' && d.type === 'BD' && d.atoms.length === 2 && d.atoms.every(t => ringSet.has(t.n))) tag = toSigma ? 'pi_sigma' : 'other';
        else if (xb.type === 'atom' && d.type === 'BD' && d.atoms.some(t => t.n === xb.Y)) tag = 'other';
        if (!tag) continue;
        contrib[tag] += e.e2; rows.push(Object.assign({ tag }, e));
      }
      const core = xb.type === 'atom' ? contrib.lp_sigma + (useRY ? contrib.lp_ry : 0) : contrib.pi_sigma;
      return { id: xb.id, label: xb.label, ionicX: !hasSigmaStar && rows.some(r => r.acceptor.type === 'LP*'), E2_core: core, E2_total: core + contrib.other + (useRY ? 0 : contrib.lp_ry), contrib, rows, verified: core >= thr, threshold: thr, includeRY: useRY,
        npaX: nbo.npa && nbo.npa[xb.X] !== undefined ? nbo.npa[xb.X] : null, npaY: xb.Y && nbo.npa ? nbo.npa[xb.Y] ?? null : null };
    });
  }
  function fragmentNPA(manifest, nbo) {
    if (!nbo.npa || !Object.keys(nbo.npa).length) return null;
    const q = { 1: 0, 2: 0 }; manifest.atoms.forEach(a => { if (nbo.npa[a.n] !== undefined) q[a.frag === 1 ? 1 : 2] += nbo.npa[a.n]; });
    return { frag1: q[1], frag2: q[2], ct: (manifest.fragCharge[1] || 0) - q[1] };
  }

  /* ------------------------------------------------------------------ */
  /*  T4: Q-Chem ALMO-EDA (EDA2) output parsing                          */
  /* ------------------------------------------------------------------ */
  function parseEDA(text) {
    const res = { units: 'kJ/mol' }; const lines = text.split(/\r?\n/);
    let iS = -1; lines.forEach((l, i) => { if (/Simplified EDA Summary/i.test(l)) iS = i; });
    const put = (k, v) => { if (!isNaN(v)) res[k] = v; };
    if (iS >= 0) {
      if (/kcal\/mol/i.test(lines[iS])) res.units = 'kcal/mol';
      for (let i = iS + 1; i < Math.min(lines.length, iS + 30); i++) {
        const m = /^\s*\[?\s*(PREPARATION|SOLVATION|FROZEN|PAULI|ELEC|DISP|POLARIZATION|CHARGE TRANSFER|TOTAL)\s+(-?\d+\.\d+)/i.exec(lines[i]);
        if (m) { const k = { PREPARATION: 'prep', SOLVATION: 'sol', FROZEN: 'frz', PAULI: 'pauli', ELEC: 'elec', DISP: 'disp', POLARIZATION: 'pol', 'CHARGE TRANSFER': 'ct', TOTAL: 'total' }[m[1].toUpperCase()]; put(k, parseFloat(m[2])); if (k === 'total') break; }
        // with ECPs Q-Chem cannot split ELEC and PAULI: "[ELEC + PAULI = 8.0364, DISP = -7.8793]"
        const mb = /\[\s*ELEC\s*\+\s*PAULI\s*=\s*(-?\d+\.\d+)\s*,\s*DISP\s*=\s*(-?\d+\.\d+)\s*\]/i.exec(lines[i]);
        if (mb) { put('elecPauli', parseFloat(mb[1])); put('disp', parseFloat(mb[2])); }
      }
    }
    const grab = (re, k) => { let v = null; lines.forEach(l => { const m = re.exec(l); if (m) { v = parseFloat(m[2]); if (/kcal/i.test(m[1])) res.units = 'kcal/mol'; } }); if (v !== null && res[k] === undefined) res[k] = v; };
    grab(/E_prp\s*\((k\w+\/mol)\)\s*=\s*(-?\d+\.\d+)/, 'prep'); grab(/E_frz\s*\((k\w+\/mol)\)\s*=\s*(-?\d+\.\d+)/, 'frz'); grab(/E_pol\s*\((k\w+\/mol)\)\s*=\s*(-?\d+\.\d+)/, 'pol');
    grab(/E_vct\s*\((k\w+\/mol)\)\s*=\s*(-?\d+\.\d+)/, 'ct'); grab(/E_ct\s*\((k\w+\/mol)\)\s*=\s*(-?\d+\.\d+)/, 'ct'); grab(/E_int\s*\((k\w+\/mol)\)\s*=\s*(-?\d+\.\d+)/, 'total');
    grab(/E_elec\s*\(ELEC\)\s*\((k\w+\/mol)\)\s*=\s*(-?\d+\.\d+)/, 'elec'); grab(/E_pauli\s*\(PAULI\)\s*\((k\w+\/mol)\)\s*=\s*(-?\d+\.\d+)/, 'pauli'); grab(/E_disp\s*\(DISP\)\s*\((k\w+\/mol)\)\s*=\s*(-?\d+\.\d+)/, 'disp');
    if (res.frz === undefined && res.pol === undefined && res.ct === undefined) return { error: 'No EDA2 energy terms found (expected "Simplified EDA Summary" or E_frz/E_pol/E_vct lines).' };
    if (res.elecPauli === undefined && res.elec !== undefined && res.pauli !== undefined) res.elecPauli = res.elec + res.pauli;
    if (res.total === undefined) res.total = (res.prep || 0) + (res.sol || 0) + (res.frz || 0) + (res.pol || 0) + (res.ct || 0);
    res.terminated = /Thank you very much for using Q-Chem/.test(text);
    return toKJ(res);
  }
  function toKJ(r) { if (r.units === 'kcal/mol') { ['prep', 'sol', 'frz', 'pauli', 'elec', 'elecPauli', 'disp', 'pol', 'ct', 'total'].forEach(k => { if (r[k] !== undefined) r[k] *= 4.184; }); r.units = 'kJ/mol'; } return r; }

  /* ------------------------------------------------------------------ */
  /*  T5: BondX Quality Index (BXQI) and mechanistic classification      */
  /* ------------------------------------------------------------------ */
  function edaDescriptors(eda, P) {
    if (!eda) return null; const e = (P ? fullParams(P) : DEFAULTS).eda;
    const pol = Math.min(0, eda.pol || 0), ct = Math.min(0, eda.ct || 0);
    const chiCT = (pol + ct) !== 0 ? ct / (pol + ct) : null;
    const attractive = [eda.elec, eda.disp, eda.pol, eda.ct].filter(v => v !== undefined && v < 0).reduce((s, v) => s + v, 0);
    const shares = attractive < 0 ? { elec: Math.min(0, eda.elec || 0) / attractive, disp: Math.min(0, eda.disp || 0) / attractive, pol: pol / attractive, ct: ct / attractive } : null;
    let mech = 'undetermined';
    if (chiCT !== null) mech = chiCT >= e.chiHigh ? 'Charge-transfer driven (n→σ* dominant)' : chiCT <= e.chiLow ? 'Polarization driven (induction dominant)' : 'Mixed polarization / charge transfer';
    let character = null;
    if (shares) { const top = Object.entries(shares).sort((a, b) => b[1] - a[1])[0]; character = { elec: 'electrostatic σ-hole', disp: 'dispersion-stabilised', pol: 'polarization', ct: 'charge-transfer' }[top[0]] + ` (${(top[1] * 100).toFixed(0)}% of attraction)`; }
    return { chiCT, ctPolRatio: pol !== 0 ? ct / pol : null, shares, mechanism: mech, character };
  }
  function computeBXQI(c, P) {
    const b = fullParams(P).bxqi; const w = b.weights; const gate = fullParams(P).nbo.threshold;
    const tiers = {};
    if (c.G != null) tiers.G = Math.max(0, Math.min(1, c.G));
    if (c.occupancy != null) tiers.D = Math.max(0, Math.min(1, c.Dscore != null ? c.Dscore : c.occupancy));
    if (c.E2 != null) tiers.E = 1 - Math.exp(-Math.max(0, c.E2) / b.E0);
    if (c.eda && c.eda.total != null) tiers.C = 1 - Math.exp(-Math.max(0, -c.eda.total) / b.C0);
    const keys = Object.keys(tiers).filter(k => w[k] > 0); if (!keys.length) return null;
    const wsum = keys.reduce((s, k) => s + w[k], 0);
    const score = 100 * Math.exp(keys.reduce((s, k) => s + (w[k] / wsum) * Math.log(Math.max(1e-4, tiers[k])), 0));
    let grade = score >= b.gradeA ? 'A · strong, high-confidence XB' : score >= b.gradeB ? 'B · moderate XB' : score >= b.gradeC ? 'C · weak XB' : 'D · contact only / doubtful';
    const electronicallyVerified = 'E' in tiers && c.E2 >= gate;
    if (b.requireGateForA) {
      if (!electronicallyVerified && score >= b.gradeA) grade = 'B · strong geometry, not electronically verified';
      if ('E' in tiers && !electronicallyVerified) grade = (score >= b.gradeC ? 'C' : 'D') + ' · rejected at NBO gate';
    }
    let consistency = null;
    if (c.E2 != null && c.eda && c.eda.ct != null) {
      const e2kJ = c.E2 * 4.184; const ct = -c.eda.ct;
      consistency = (e2kJ > b.consE2hi && ct < b.consCTlo) ? 'Inconsistent: large NBO E(2) but negligible ALMO CT — check basis-set superposition (enable EDA_BSSE) or NBO delocalisation artefact.'
        : (e2kJ < b.consE2lo && ct > b.consCThi) ? 'Inconsistent: ALMO CT large but NBO E(2) small — CT likely into Rydberg/diffuse space; inspect COVPs.' : 'Consistent orbital (NBO) and variational (ALMO) charge-transfer pictures.';
    }
    return { score, grade, tiers, electronicallyVerified, provenance: ['G', 'D', 'E', 'C'].filter(k => k in tiers).join('·'), consistency, weights: w };
  }

  /* ------------------------------------------------------------------ */
  /*  Utilities                                                          */
  /* ------------------------------------------------------------------ */
  function toCSV(rows, cols) {
    const esc = v => { if (v === null || v === undefined) return ''; const s = typeof v === 'number' ? (Number.isInteger(v) ? String(v) : v.toFixed(4)) : String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
    return [cols.map(c => c.h).join(','), ...rows.map(r => cols.map(c => esc(typeof c.k === 'function' ? c.k(r) : r[c.k])).join(','))].join('\n') + '\n';
  }

  return {
    VERSION, Z, COV, DEFAULTS, PARAM_SPEC, FF: FFDATA, defaults, mergeParams, fullParams, getPath, setPath, assumptions,
    parsePDB, structureReport, filterChains, writePDB, coordOf, atomLabel, resLabel, resId,
    detectXB, prepareDetection, detectFrame, decorate, geometricScore, parseXYZFrames, parseDCD, dropVsites, analyzeTrajectory, trajAccumulator, Dscore,
    parseCIF, parseXYZ, detectMaterialXB, buildMaterialQM,
    ffVariant, assignForceField, titratableResidues, parseMOL2, mol2ForStructure,
    ligandProtonation, formulCheck, ORCA_NBO_CAUTION, QCHEM_NBO_CAUTION, QCHEM_NBO_G_CAUTION, G16_NBO7_CAUTION, buildQMModel, embeddingCharges, writeGaussian, writeORCA, writeQChemNBO, writeQChemEDA, writeXYZ, chargesCSV,
    parseNBO, nboGate, verifyGeometry, fragmentNPA, parseEDA, edaDescriptors, computeBXQI,
    hashString, toCSV, angle, dist
  };
});
