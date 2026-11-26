#!/usr/bin/env python3
"""09_performance.py
Run time of the BondX engine (Node.js, same code as the browser app).
  1. T1 parse + detection time for every protein test structure and for tiled copies of 1P5E (up to ~74 000 atoms).
  2. T2 trajectory throughput (frames per second) for 1US0+IDD594 and 1P5E with 200 synthetic frames.
  3. T3a QM-model build time (cluster and QM/MM) from results/protein_models.csv runs.
Each timing is the median of 5 repeats. Output: results/performance.csv
"""
import os, sys, csv, json, subprocess, statistics
from bxutil import cli, HERE

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
PDB = os.path.join(ROOT, 'data', 'protein')
CHAINS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'


def tile(src, n, out):
    """n copies of a PDB, shifted by 150 A along x, with new chain letters."""
    L = [l for l in open(src) if l.startswith(('ATOM', 'HETATM'))]
    chains = sorted({l[21] for l in L}); k = 0
    with open(out, 'w') as f:
        for c in range(n):
            for l in L:
                ch = CHAINS[(c * len(chains) + chains.index(l[21])) % len(CHAINS)]
                x = float(l[30:38]) + 150.0 * c
                f.write(l[:21] + ch + l[22:30] + '%8.3f' % x + l[38:])
            k += len(L)
        f.write('END\n')
    return k


def time_detect(path, rep=5):
    t = []
    for _ in range(rep):
        d = cli('detect-pdb', path); t.append(d['msParse'] + d['msDetect']); n = d['atoms']; nx = len([h for h in d['hits'] if h['status'] == 'XB'])
    return n, statistics.median(t), nx


TRAJ_JS = r'''
const B = require(process.argv[2]); const fs = require('fs');
const S = B.parsePDB(fs.readFileSync(process.argv[3], 'utf8')); const P = B.defaults();
const f0 = S.frames[0]; let seed = 7; const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647 - 0.5; };
const frames = []; for (let k = 0; k < 200; k++) { const f = Float32Array.from(f0); for (let i = 0; i < f.length; i++) f[i] += rnd() * 0.4; frames.push(f); }
const t = []; for (let r = 0; r < 5; r++) { const acc = B.trajAccumulator(S, P, {}); const t0 = process.hrtime.bigint(); frames.forEach((f, i) => acc.add(f, null, i)); const res = acc.finish(); t.push(Number(process.hrtime.bigint() - t0) / 1e6); }
t.sort((a, b) => a - b); console.log(JSON.stringify({ atoms: S.atoms.length, frames: 200, ms: t[2] }));
'''


def main():
    rows = []
    for f in sorted(os.listdir(PDB)):
        if f.endswith('.pdb'):
            n, ms, nx = time_detect(os.path.join(PDB, f)); rows.append(dict(test='T1 detect', system=f[:-4], atoms=n, ms=round(ms, 1), xb=nx))
    for c in (2, 4, 8):
        out = '/tmp/1p5e_x%d.pdb' % c; tile(os.path.join(PDB, '1p5e.pdb'), c, out)
        n, ms, nx = time_detect(out); rows.append(dict(test='T1 detect', system='1p5e x%d' % c, atoms=n, ms=round(ms, 1), xb=nx))
    js = '/tmp/bx_traj_bench.js'; open(js, 'w').write(TRAJ_JS)
    for s in ('1us0_idd594', '1p5e'):
        r = json.loads(subprocess.run(['node', js, os.path.join(HERE, 'engine.js'), os.path.join(PDB, s + '.pdb')], capture_output=True, text=True).stdout)
        rows.append(dict(test='T2 trajectory', system=s, atoms=r['atoms'], ms=round(r['ms'], 1), xb='', frames_per_s=round(r['frames'] / r['ms'] * 1000, 1)))
    for s, lig in (('1us0_idd594', 'LDT:A:320'), ('1p5e', 'TBS:A:301'), ('4gql', 'R47:A:306')):
        for env in ('cluster', 'qmmm'):
            t = [cli('build-pdb', os.path.join(PDB, s + '.pdb'), '--ligand', lig, '--out', '/tmp/bx_perf', '--name', 'p', '--programs', 'none', '--noEda', '1', '--P', json.dumps({'qm': {'env': env}}))['ms'] for _ in range(5)]
            rows.append(dict(test='T3a build ' + env, system=s, atoms='', ms=statistics.median(t), xb=''))
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(ROOT, 'results', 'performance.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    for r in rows: print(r)


if __name__ == '__main__':
    main()
