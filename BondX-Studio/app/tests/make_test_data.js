// Builds a synthetic but chemically sensible complex for testing BondX.
const fs = require('fs');
const sub = (a, b) => a.map((v, i) => v - b[i]), add = (a, b) => a.map((v, i) => v + b[i]), scl = (a, s) => a.map(v => v * s);
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2], cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const unit = a => scl(a, 1 / Math.hypot(...a)); const rad = d => d * Math.PI / 180;
function nerf(A, B, C, bond, ang, tor) {
  const bc = unit(sub(C, B)); const n = unit(cross(sub(B, A), bc)); const m = cross(n, bc);
  const d2 = [-bond * Math.cos(rad(ang)), bond * Math.sin(rad(ang)) * Math.cos(rad(tor)), bond * Math.sin(rad(ang)) * Math.sin(rad(tor))];
  return add(C, add(add(scl(bc, d2[0]), scl(m, d2[1])), scl(n, d2[2])));
}
const seq = ['ALA', 'PHE', 'ASP', 'SER', 'MET', 'GLY'];
const phi = -65, psi = -40; const atoms = [];
let N = [0, 0, 0], CA = [1.458, 0, 0], C = nerf([0, 1, 0], N, CA, 1.525, 111.2, -60);
for (let r = 0; r < seq.length; r++) {
  const rn = seq[r]; const res = r + 10; const A = (name, p, el) => atoms.push({ name, rn, res, p, el: el || name[0], chain: 'A' });
  A('N', N); A('CA', CA); A('C', C);
  const Nn = nerf(N, CA, C, 1.329, 116.2, psi); const O = nerf(Nn, CA, C, 1.231, 120.5, 180); A('O', O);
  if (rn !== 'GLY') {
    const CB = nerf(C, N, CA, 1.53, 110.5, -122.6); A('CB', CB);
    if (rn === 'PHE') { const CG = nerf(N, CA, CB, 1.50, 113.8, -60); A('CG', CG); const CD1 = nerf(CA, CB, CG, 1.39, 120.7, 90), CD2 = nerf(CA, CB, CG, 1.39, 120.7, -90);
      A('CD1', CD1); A('CD2', CD2); const CE1 = nerf(CB, CG, CD1, 1.39, 120, 180), CE2 = nerf(CB, CG, CD2, 1.39, 120, 180); A('CE1', CE1); A('CE2', CE2); A('CZ', nerf(CG, CD1, CE1, 1.39, 120, 0)); }
    if (rn === 'ASP') { const CG = nerf(N, CA, CB, 1.52, 113, -60); A('CG', CG); A('OD1', nerf(CA, CB, CG, 1.25, 118.4, 0), 'O'); A('OD2', nerf(CA, CB, CG, 1.25, 118.4, 180), 'O'); }
    if (rn === 'SER') A('OG', nerf(N, CA, CB, 1.42, 111, 60), 'O');
    if (rn === 'MET') { const CG = nerf(N, CA, CB, 1.52, 114, -60); A('CG', CG); const SD = nerf(CA, CB, CG, 1.81, 112.7, 180); A('SD', SD, 'S'); A('CE', nerf(CB, CG, SD, 1.79, 100.8, 70)); }
  }
  if (r === seq.length - 1) A('OXT', nerf(Nn, CA, C, 1.25, 117, 0), 'O');
  const CAn = nerf(CA, C, Nn, 1.458, 121.7, 180); const Cn = nerf(C, Nn, CAn, 1.525, 111.2, phi);
  N = Nn; CA = CAn; C = Cn;
}
const get = (res, name) => atoms.find(a => a.res === res && a.name === name).p;
function benzene(X, dir, el, xc, rn, res, perp) {
  const C1 = add(X, scl(dir, xc)); const cen = add(C1, scl(dir, 1.39));
  const p = unit(sub(perp, scl(dir, dot(perp, dir))));
  const out = [{ name: el.length === 2 ? el + '1' : el + '1', rn, res, p: X, el, chain: 'L', het: true }];
  for (let k = 0; k < 6; k++) { const th = rad(60 * k); out.push({ name: 'C' + (k + 1), rn, res, p: add(cen, add(scl(dir, -1.39 * Math.cos(th)), scl(p, 1.39 * Math.sin(th)))), el: 'C', chain: 'L', het: true }); }
  return out;
}
// Iodobenzene halogen-bonded to ALA10 backbone O (d = 3.00 Å, C-I...O = 178°, I...O=C ~ 125°)
const O = get(10, 'O'), Cc = get(10, 'C'); const co = unit(sub(O, Cc)); const ax = unit(cross(co, [0.3, 0.2, 1]));
const rot = (v, k, a) => { const c = Math.cos(a), s = Math.sin(a); return add(add(scl(v, c), scl(cross(k, v), s)), scl(k, dot(k, v) * (1 - c))); };
let best = null;
for (let tw = 0; tw < 360; tw += 10) { // choose a clash-free approach direction ~125° from C=O
  const d0 = rot(rot(co, ax, rad(55)), co, rad(tw)); const Ii = add(O, scl(d0, 3.0));
  const lig = benzene(Ii, d0, 'I', 2.10, 'X', 0, ax); const mind = Math.min(...lig.flatMap(l => atoms.map(a => Math.hypot(...sub(a.p, l.p)))));
  if (!best || mind > best.m) best = { d0, m: mind };
}
const dI = best.d0; const I = add(O, scl(dI, 3.00));
atoms.push(...benzene(I, dI, 'I', 2.10, 'IBZ', 901, ax));
// Bromobenzene over PHE11 ring (Br...centroid 3.45 Å along the ring normal)
const ring = ['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2'].map(n => get(11, n)); const cen = scl(ring.reduce(add), 1 / 6);
let nrm = unit(cross(sub(ring[1], ring[0]), sub(ring[2], ring[0]))); if (dot(nrm, sub(cen, get(11, 'CB'))) < 0) { }
const Br = add(cen, scl(nrm, 3.45)); atoms.push(...benzene(Br, nrm, 'BR', 1.90, 'BBZ', 902, [1, 0, 0]));
// waters
atoms.push({ name: 'O', rn: 'HOH', res: 1001, p: add(get(13, 'OG'), [2.8, 0, 0]), el: 'O', chain: 'W', het: true });
function pdb(list, title) {
  const L = ['HEADER    TEST COMPLEX                            18-NOV-26   TEST', 'TITLE     BONDX SYNTHETIC XB TEST', 'COMPND    MOL_ID: 1;', 'COMPND   2 MOLECULE: TEST PEPTIDE;', 'COMPND   3 CHAIN: A;'];
  list.forEach((a, i) => { const nm = a.el.length === 1 && a.name.length < 4 ? ' ' + a.name.padEnd(3) : a.name.padEnd(4);
    L.push((a.het ? 'HETATM' : 'ATOM  ') + String(i + 1).padStart(5) + ' ' + nm + ' ' + a.rn.padStart(3) + ' ' + a.chain + String(a.res).padStart(4) + '    ' + a.p.map(v => v.toFixed(3).padStart(8)).join('') + '  1.00 20.00          ' + a.el.padStart(2)); });
  L.push('END'); return L.join('\n') + '\n';
}
fs.mkdirSync(__dirname + '/data', { recursive: true });
fs.writeFileSync(__dirname + '/data/demo_complex.pdb', pdb(atoms));
// trajectory: 60 frames, thermal jitter + slow drift of the Br ligand away after frame 30
let traj = ''; let seed = 7; const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647 - 0.5; };
const frames = [];
for (let f = 0; f < 60; f++) {
  const fr = atoms.map(a => { let p = add(a.p, [rnd() * 0.25, rnd() * 0.25, rnd() * 0.25]); if (a.rn === 'BBZ' && f > 30) p = add(p, scl(nrm, 0.08 * (f - 30))); return Object.assign({}, a, { p }); });
  frames.push(fr); traj += 'MODEL     ' + String(f + 1).padStart(4) + '\n' + pdb(fr).split('\n').filter(l => /^(ATOM|HETATM)/.test(l)).join('\n') + '\nENDMDL\n';
}
fs.writeFileSync(__dirname + '/data/demo_traj.pdb', traj + 'END\n');
// DCD (CHARMM format, little endian, no unit cell)
const n = atoms.length; const parts = [];
const rec = (buf) => { const h = Buffer.alloc(4); h.writeInt32LE(buf.length); parts.push(h, buf, h); };
const hdr = Buffer.alloc(84); hdr.write('CORD', 0); hdr.writeInt32LE(frames.length, 4); hdr.writeInt32LE(0, 8); hdr.writeInt32LE(1, 12); hdr.writeInt32LE(24, 80); rec(hdr);
const tb = Buffer.alloc(4 + 80); tb.writeInt32LE(1, 0); tb.write('BONDX TEST DCD', 4); rec(tb);
const nb = Buffer.alloc(4); nb.writeInt32LE(n); rec(nb);
for (const fr of frames) for (let c = 0; c < 3; c++) { const b = Buffer.alloc(4 * n); fr.forEach((a, k) => b.writeFloatLE(a.p[c], 4 * k)); rec(b); }
fs.writeFileSync(__dirname + '/data/demo_traj.dcd', Buffer.concat(parts));
console.log('atoms', atoms.length, 'I...O', Math.hypot(...sub(I, O)).toFixed(3), 'Br...cen', Math.hypot(...sub(Br, cen)).toFixed(3));
