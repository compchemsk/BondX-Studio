#!/usr/bin/env python3
"""15_hiyd_analysis.py
BondX re-analysis of the hIYD I-Tyr unbinding data (Karmakar & Mishra, J. Chem. Inf. Model. 2026, 66, 11191).
No new MD or QM calculation is done. Inputs (hIYD_for_BondX):
  US/C1..C3/traj/*.pdb,*.dcd   stripped umbrella-sampling windows (36 per pathway, water removed)
  US/C1..C3/halo.out           published frame-wise structural screen (first 10 ns = 500 frames per window)
  US/C1..C3/qm-halo/*-nbo.log  published QM-cluster NBO outputs (Gaussian 16), one representative pose per residue
  US/C1..C3/qm-halo/role-halo.dat  window, frame, d, theta and the published NBO value per residue
  collect.log                  frame stride of every stripped trajectory
Steps
  1. node hiyd_bondx.js  -> results/hiyd/bondx_frames.jsonl (T1 contacts of every frame)
  2. node hiyd_nbo.js    -> results/hiyd/nbo_clusters.json (BondX NBO parser and gate on the 56 cluster outputs)
  3. exact check: BondX T1 contacts vs the published screen on the frames present in both
  4. occupancy per window and per residue (all 25 ns); acceptor classes; strict BondX filter
  5. NBO re-analysis with a second, independent parser (sections 'within unit' / 'from unit')
Usage: python3 15_hiyd_analysis.py <project root> <hIYD_for_BondX root>
Pathway C1, C2, C3 = Pathway A, B, C of the paper. Residues: MD numbers 1-440 (1-220 chain A = PDB 71-290; 221-440 chain B = PDB 71*-290*), 441/443 FRA/FRB, 442/444 MIA/MIB.
Chain A PDB = n + 70, chain B PDB = n - 150 (*).
The published screen (MDTraj on the parm7) writes 0-based residue indices: ARG93 there is Amber residue 94 = R164.
"""
import os, sys, re, json, csv, math, subprocess
from collections import defaultdict, Counter
import numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
DATA = sys.argv[2] if len(sys.argv) > 2 else 'hIYD_for_BondX'  # folder with the hIYD trajectories (give it as the second argument)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, 'results', 'hiyd'); os.makedirs(OUT, exist_ok=True)
PATH = {'C1': 'A', 'C2': 'B', 'C3': 'C'}
WINS = [round(3.8 + 0.2 * i, 1) for i in range(26)] + [9.8 + i for i in range(10)]
WINS = [round(w, 1) for w in WINS]
CUT = {'N': 3.97, 'O': 3.82, 'S': 4.15, 'AR': 4.00}

def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode: sys.exit(r.stderr)
    return r.stdout.strip()

def pdbnum(res, n):
    if res in ('FRA', 'FRB'): return 'FMN'
    if res in ('MIA', 'MIB'): return 'I-Tyr(2)'
    return f'{n + 70}' if n <= 220 else f'{n - 150}*'

ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'HIE': 'H', 'HID': 'H', 'HIP': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
def label(res, n):
    if res in ('FRA', 'FRB'): return 'FMN'
    return ONE.get(res, res) + pdbnum(res, n)
def amber_name(res):  # halo.out writes HIE/HID as HIS
    return 'HIS' if res in ('HIE', 'HID', 'HIP') else res

# BondX Lewis classes (engine LEWIS table); everything else 'strong'
LEWIS_WEAK = {('_BB', 'N'), ('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('LYS', 'NZ'), ('TRP', 'NE1'), ('ASN', 'ND2'), ('GLN', 'NE2'), ('PRO', 'N')}
def acc_class(res, atom, typ):
    if typ in ('aromatic', 'pi'): return 'aromatic ring'
    if res in ('FRA', 'FRB', 'MIB'): return 'FMN ' + atom[0]
    if atom == 'N': return 'backbone N-H'
    if atom in ('O', 'OXT'): return 'backbone C=O'
    if (res, atom) in LEWIS_WEAK or (res == 'HIS' and False): return 'side-chain N-H (no lone pair)'
    if atom[0] == 'N': return 'side-chain N (lone pair)'
    if atom[0] == 'S': return 'thioether S'
    return 'side-chain O'
def lewis(res, atom, typ):
    if typ in ('aromatic', 'pi'): return 'strong'
    if atom == 'N' and res not in ('FRA', 'FRB', 'MIB'): return 'weak'
    return 'weak' if (amber_name(res), atom) in LEWIS_WEAK else 'strong'

# ------------------------------------------------------------------ 1-2. BondX runs
if not (len(sys.argv) > 3 and sys.argv[3] == '--no-node'):
    print(run(['node', os.path.join(HERE, 'hiyd_bondx.js'), DATA, os.path.join(OUT, 'bondx_frames.jsonl')]))
    print(run(['node', os.path.join(HERE, 'hiyd_nbo.js'), DATA, os.path.join(OUT, 'nbo_clusters.json')]))

# ------------------------------------------------------------------ strides from collect.log
stride, nraw = {}, {}
for l in open(os.path.join(DATA, 'collect.log')):
    m = re.match(r'\s+(C\d)_US-([\d.]+): (\d+) frames, stride (\d+)', l)
    if m: stride[(m.group(1), float(m.group(2)))] = int(m.group(4)); nraw[(m.group(1), float(m.group(2)))] = int(m.group(3))
assert len(stride) == 108, len(stride)

# ------------------------------------------------------------------ BondX frames
frames = defaultdict(dict)   # (C,w) -> k -> hits
nfr = {}
for l in open(os.path.join(OUT, 'bondx_frames.jsonl')):
    r = json.loads(l); frames[(r['C'], r['w'])][r['k']] = r['hits']; nfr[(r['C'], r['w'])] = r['n']
assert len(frames) == 108
for key, n in nfr.items():
    assert n == math.ceil(nraw[key] / stride[key]) and len(frames[key]) == n, (key, n)

# ------------------------------------------------------------------ published screen (halo.out)
pub = defaultdict(lambda: defaultdict(list))  # (C,w) -> frame -> [(res,num,type,atom,d_A,theta)]
LINE = re.compile(r'^\s+([A-Z0-9+\-]{2,4}?)(\d+) \((N|O|S|aromatic), ([^)]*)\), Distance: ([\d.]+) nm, Donor Pos: \[[^\]]*\], Angle: ([\d.]+)')
for C in PATH:
    w = fr = None
    txt = open(os.path.join(DATA, 'US', C, 'halo.out')).read()
    txt = re.sub(r'\[[^\]]*\]', '[]', txt, flags=re.S)   # numpy may wrap the position over two lines
    for l in txt.split('\n'):
        m = re.match(r'For trajectory ([\d.]+)', l)
        if m: w = float(m.group(1)); fr = None; continue
        m = re.match(r'Frame (\d+):', l)
        if m: fr = int(m.group(1)); continue
        m = LINE.match(l)
        if m:
            res, num, typ, atom = m.group(1), int(m.group(2)), m.group(3), m.group(4).strip()
            pub[(C, w)][fr].append((res, num + 1, typ, 'ring' if typ == 'aromatic' else atom, float(m.group(5)) * 10, float(m.group(6))))
npub_lines = sum(len(v) for d in pub.values() for v in d.values())

# ------------------------------------------------------------------ 3. exact frame check
def bx_key(h):
    res = amber_name(h['res'])
    return (res, h['seq'], 'aromatic' if h['type'] == 'pi' else h['el'], 'ring' if h['type'] == 'pi' else h['atom'])
chk = Counter(); dmax = tmax = 0.0; missing, extra = [], []
for (C, w), byk in frames.items():
    s = stride[(C, w)]
    for k, hits in byk.items():
        f = k * s
        if f > 499: continue
        chk['frames'] += 1
        P = [p for p in pub[(C, w)].get(f, []) if p[0] != 'HOH']
        chk['water_lines'] += len(pub[(C, w)].get(f, [])) - len(P)
        Bk = defaultdict(list)
        for h in hits: Bk[bx_key(h)].append(h)
        used = set()
        for p in P:
            key = (p[0], p[1], p[2], p[3]); chk['pub_lines'] += 1
            cand = [h for h in Bk.get(key, []) if id(h) not in used]
            if not cand: missing.append((C, w, f, p)); continue
            h = min(cand, key=lambda h: abs(h['d'] - p[4]) + abs(h['t1'] - p[5]) / 100)
            used.add(id(h)); chk['matched'] += 1
            dmax = max(dmax, abs(h['d'] - p[4])); tmax = max(tmax, abs(h['t1'] - p[5]))
        for key, hs in Bk.items():
            for h in hs:
                if id(h) not in used: extra.append((C, w, f, key, h['d'], h['t1']))
chk.update({'missing': len(missing), 'extra': len(extra)})
check = {'framesCompared': chk['frames'], 'publishedLinesNoWater': chk['pub_lines'], 'waterLinesSkipped': chk['water_lines'], 'matched': chk['matched'],
         'missing': chk['missing'], 'extra': chk['extra'], 'maxAbsDeltaD_A': round(dmax, 4), 'maxAbsDeltaTheta_deg': round(tmax, 3),
         'missingList': [list(map(str, m)) for m in missing[:20]], 'extraList': [list(map(str, e)) for e in extra[:20]]}
print('frame check', json.dumps({k: v for k, v in check.items() if not k.endswith('List')}))

# ------------------------------------------------------------------ 4. occupancy per window
def strict(h):  # BondX strict: Lewis-base acceptor; ring contacts within 60 deg of the ring normal
    lw = lewis(h['res'], h['atom'], h['type'])
    if h['type'] == 'pi': return h['t2'] is not None and h['t2'] <= 60
    return lw == 'strong'
rows = []
for C in PATH:
    for w in WINS:
        byk = frames[(C, w)]; n = len(byk); s = stride[(C, w)]
        anyA = sum(1 for hs in byk.values() if hs); anyS = sum(1 for hs in byk.values() if any(strict(h) for h in hs))
        k10 = [k for k in byk if k * s <= 499]
        any10 = sum(1 for k in k10 if byk[k])
        P = pub[(C, w)]
        pubAll = sum(1 for f, v in P.items() if v); pubNoW = sum(1 for f, v in P.items() if any(p[0] != 'HOH' for p in v))
        rows.append({'pathway': PATH[C], 'C': C, 'xi': w, 'stride': s, 'frames_25ns': n, 'bondx_T1_frac': anyA / n, 'bondx_strict_frac': anyS / n,
                     'frames_10ns_in_dcd': len(k10), 'bondx_T1_frac_10ns': any10 / len(k10), 'pub_frac_10ns': pubAll / 500, 'pub_frac_10ns_nowater': pubNoW / 500})
with open(os.path.join(OUT, 'hiyd_window_occupancy.csv'), 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)

# ------------------------------------------------------------------ acceptor classes and residues (all frames)
cls_count = {C: Counter() for C in PATH}; cls_strict = {C: Counter() for C in PATH}
res_frames = {C: defaultdict(set) for C in PATH}; res_frames_strict = {C: defaultdict(set) for C in PATH}
nframes_C = {C: 0 for C in PATH}; xi_by_res = {C: defaultdict(set) for C in PATH}
G_all, G_strict = [], []
for (C, w), byk in frames.items():
    nframes_C[C] += len(byk)
    for k, hs in byk.items():
        for h in hs:
            typ = 'pi' if h['type'] == 'pi' else h['el']
            c = acc_class(amber_name(h['res']), h['atom'], typ); cls_count[C][c] += 1
            lab = label(h['res'], h['seq']); res_frames[C][lab].add((w, k)); xi_by_res[C][lab].add(w); G_all.append(h['G'])
            if strict(h): cls_strict[C][c] += 1; res_frames_strict[C][lab].add((w, k)); G_strict.append(h['G'])
classes = sorted({c for C in PATH for c in cls_count[C]}, key=lambda c: (-sum(cls_count[C][c] for C in PATH), c))
with open(os.path.join(OUT, 'hiyd_acceptor_classes.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['class'] + [f'contacts_{PATH[C]}' for C in PATH] + ['total', 'kept_by_strict_filter'])
    for c in classes: wr.writerow([c] + [cls_count[C][c] for C in PATH] + [sum(cls_count[C][c] for C in PATH), sum(cls_strict[C][c] for C in PATH)])
res_rows = []
for C in PATH:
    for lab, fs in res_frames[C].items():
        res_rows.append({'pathway': PATH[C], 'residue': lab, 'frames_T1': len(fs), 'frames_strict': len(res_frames_strict[C].get(lab, ())),
                         'windows': len({w for w, k in fs}), 'xi_min': min(w for w, k in fs), 'xi_max': max(w for w, k in fs)})
res_rows.sort(key=lambda r: (r['pathway'], -r['frames_T1']))
with open(os.path.join(OUT, 'hiyd_residue_contacts.csv'), 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(res_rows[0])); wr.writeheader(); wr.writerows(res_rows)

# ------------------------------------------------------------------ 5. NBO re-analysis
NBL = re.compile(r'^\s*(\d+|\*+)\.\s*(BD\*?|LP\*?|CR|RY\*?|LV|3C\*?)\s*\(\s*(\d+)\s*\)\s*([A-Z][a-z]?)\s*(\d+)(?:\s*-\s*([A-Z][a-z]?)\s*(\d+))?(?:\s*-\s*([A-Z][a-z]?)\s*(\d+))?\s*/\s*(\d+|\*+)\.\s*(BD\*?|LP\*?|CR|RY\*?|LV|3C\*?)\s*\(\s*(\d+)\s*\)\s*([A-Z][a-z]?)\s*(\d+)(?:\s*-\s*([A-Z][a-z]?)\s*(\d+))?(?:\s*-\s*([A-Z][a-z]?)\s*(\d+))?\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')
def py_nbo(fn, iI, iC):
    """Independent parser: LP -> sigma*(C-I) entries with their NBO section ('within unit u' or 'from unit u to unit v')."""
    L = open(fn).read().split('\n')
    st = max(i for i, l in enumerate(L) if 'SECOND ORDER PERTURBATION THEORY ANALYSIS' in l.upper())
    sec, out = None, []
    for l in L[st:]:
        if 'NATURAL BOND ORBITALS (Summary)' in l: break
        m = re.match(r'\s*within unit\s+(\d+)', l)
        if m: sec = ('within', int(m.group(1))); continue
        m = re.match(r'\s*from unit\s+(\d+) to unit\s+(\d+)', l)
        if m: sec = ('from', int(m.group(1)), int(m.group(2))); continue
        m = NBL.match(l)
        if not m: continue
        g = m.groups()
        acc_atoms = {int(g[13])} | ({int(g[15])} if g[15] else set())
        if g[10] == 'BD*' and acc_atoms == {iI, iC} and g[1] == 'LP':
            out.append({'donorEl': g[3], 'donor': int(g[4]), 'e2': float(g[18]), 'sec': sec})
    return out

clusters = json.load(open(os.path.join(OUT, 'nbo_clusters.json')))
role = {}
for C in PATH:
    for l in open(os.path.join(DATA, 'US', C, 'qm-halo', 'role-halo.dat')):
        p = l.split()
        if p: role[(C, p[0])] = {'xi': float(p[1]), 'frame': int(p[2]), 'd_nm': float(p[3]), 'theta': float(p[4]), 'E_pub': float(p[7])}
nbo_rows = []
for r in clusters:
    C, resid = r['C'], r['residue']; ro = role[(C, resid)]
    m = re.match(r'([A-Z]{3})(\d+)', resid); rn, num = m.group(1), int(m.group(2)) + 1   # 0-based index in the published files
    # structural acceptor from the published screen: same window, frame, residue, distance and angle
    hits = [p for p in pub[(C, ro['xi'])].get(ro['frame'], []) if p[0] == rn and p[1] == num and abs(p[4] - ro['d_nm'] * 10) < 0.006 and abs(p[5] - ro['theta']) < 0.006]
    assert len(hits) >= 1, (C, resid, ro)
    p = hits[0]; typ = p[2]; atom = p[3]
    # cluster T1 pair(s)
    pairs = [(g, ga) for g, ga in zip(r['geo'], r['gate']) if g['dok'] and g['t1ok']]
    want = 'AR' if typ == 'aromatic' else typ
    pm = [(g, ga) for g, ga in pairs if g['el'] == want and abs(g['d'] - p[4]) < 0.02 and abs(g['t1'] - p[5]) < 0.2]
    assert len(pm) == 1, (C, resid, typ, p, [(g['el'], g['d'], g['t1']) for g, ga in pairs])
    g, ga = pm[0]
    # second parser
    ent = py_nbo(os.path.join(DATA, 'US', C, 'qm-halo', resid + '-nbo.log'), r['I'], r['Cx'])
    unitI = None
    allsec = [e['sec'] for e in ent]
    # the unit that holds the iodine is the target of every 'from' section and the unit of 'within' entries of fragment-1 donors
    intra_py = sum(e['e2'] for e in ent if e['sec'][0] == 'within')
    inter_py = sum(e['e2'] for e in ent if e['sec'][0] == 'from')
    fr1 = {x['donor']: x['frag'] for x in r['lpRows']}
    node_intra = sum(x['e2'] for x in r['lpRows'] if x['frag'] == 1); node_inter = sum(x['e2'] for x in r['lpRows'] if x['frag'] == 2)
    assert abs(intra_py - node_intra) < 1e-9 and abs(inter_py - node_inter) < 1e-9, (C, resid, intra_py, node_intra, inter_py, node_inter)
    ons = [e for e in ent if e['donorEl'] in ('O', 'N', 'S')]
    emax = max([e['e2'] for e in ons], default=0.0)
    src = sorted({('intra' if e['sec'][0] == 'within' else 'inter') for e in ons if e['e2'] == emax}) or ['none']
    lw = lewis(rn, atom, typ)
    nbo_rows.append({'pathway': PATH[C], 'residue_amber': resid, 'residue': label(rn, num), 'xi': ro['xi'], 'frame': ro['frame'],
                     'acceptor_atom': atom, 'acceptor_type': typ, 'acceptor_class': acc_class(rn, atom, typ), 'lewis': lw,
                     'd_A': round(p[4], 2), 'theta1': round(p[5], 2), 'E_pub': ro['E_pub'], 'E_max_line': round(emax, 2), 'E_pub_source': '+'.join(src),
                     'E2_intra_LP': round(node_intra, 2), 'E2_intra_phenolateO': round(sum(x['e2'] for x in r['lpRows'] if x['frag'] == 1 and x['donor'][0] == 'O'), 2), 'E2_inter_LP': round(node_inter, 2), 'E2_bondx_pair': round(ga['E2_core'], 2), 'E2_RY_pair': round(ga.get('E2_RY') or 0, 2), 'hb_LPI_to_YH': round(ga.get('hbY') or 0, 2), 'LP_on_acceptor': ga.get('lpOnY'), 'verdict_v3': ga.get('verdict'),
                     'bondx_verified': ga['E2_core'] >= 0.5, 'bondx_verified_lewis': ga['E2_core'] >= 0.5 and lw == 'strong', 'nAtoms': r['nAtoms']})
with open(os.path.join(OUT, 'hiyd_nbo_reanalysis.csv'), 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(nbo_rows[0])); wr.writeheader(); wr.writerows(nbo_rows)

# ------------------------------------------------------------------ summary
def q(x): return round(float(x), 3)
S = {'windows': len(frames), 'frames_total': sum(len(v) for v in frames.values()), 'frames_per_window_min': min(nfr.values()), 'frames_per_window_max': max(nfr.values()),
     'stride_values': sorted(set(stride.values())), 'pub_lines_total': npub_lines, 'check': {k: v for k, v in check.items() if not k.endswith('List')}}
for C in PATH:
    R = [r for r in rows if r['C'] == C]
    S[f'occ_{PATH[C]}'] = {'windows_with_T1': sum(1 for r in R if r['bondx_T1_frac'] > 0), 'windows_with_strict': sum(1 for r in R if r['bondx_strict_frac'] > 0),
        'mean_T1_xi_le_8_8': q(np.mean([r['bondx_T1_frac'] for r in R if r['xi'] <= 8.8])), 'mean_T1_xi_ge_9_8': q(np.mean([r['bondx_T1_frac'] for r in R if r['xi'] >= 9.8])),
        'mean_strict_xi_le_8_8': q(np.mean([r['bondx_strict_frac'] for r in R if r['xi'] <= 8.8])), 'mean_strict_xi_ge_9_8': q(np.mean([r['bondx_strict_frac'] for r in R if r['xi'] >= 9.8])),
        'max_T1': q(max(r['bondx_T1_frac'] for r in R)), 'xi_max_T1': max(R, key=lambda r: r['bondx_T1_frac'])['xi'],
        'max_strict': q(max(r['bondx_strict_frac'] for r in R)), 'xi_max_strict': max(R, key=lambda r: r['bondx_strict_frac'])['xi'],
        'first_xi_strict_ge_0_1': next((r['xi'] for r in R if r['bondx_strict_frac'] >= 0.1), None),
        'contacts_T1': sum(cls_count[C].values()), 'contacts_strict': sum(cls_strict[C].values()), 'frames': nframes_C[C],
        'residues_T1': len(res_frames[C]), 'residues_strict': len(res_frames_strict[C]),
        'top_strict': sorted(((lab, len(v)) for lab, v in res_frames_strict[C].items()), key=lambda x: -x[1])[:6]}
tot = Counter(); tots = Counter()
for C in PATH: tot.update(cls_count[C]); tots.update(cls_strict[C])
S['classes'] = {c: [tot[c], tots[c]] for c in classes}
S['contacts_total'] = sum(tot.values()); S['contacts_strict_total'] = sum(tots.values())
S['G_median_all'] = q(np.median(G_all)); S['G_median_strict'] = q(np.median(G_strict))
weakN = sum(tot[c] for c in classes if c in ('backbone N-H', 'side-chain N-H (no lone pair)'))
S['frac_contacts_no_lone_pair'] = q(weakN / sum(tot.values()))
N = len(nbo_rows)
S['nbo'] = {'clusters': N, 'per_pathway': dict(Counter(r['pathway'] for r in nbo_rows)), 'all_terminated': all(r['terminated'] for r in clusters),
    'pub_equals_max_line': sum(1 for r in nbo_rows if abs(r['E_pub'] - r['E_max_line']) < 0.006),
    'pub_source_intra': sum(1 for r in nbo_rows if r['E_pub_source'] == 'intra'), 'pub_source_inter': sum(1 for r in nbo_rows if r['E_pub_source'] == 'inter'),
    'pub_source_none': sum(1 for r in nbo_rows if r['E_pub_source'] == 'none'),
    'pub_ge_0_5': sum(1 for r in nbo_rows if r['E_pub'] >= 0.5),
    'bondx_verified': sum(1 for r in nbo_rows if r['bondx_verified']), 'bondx_verified_lewis': sum(1 for r in nbo_rows if r['bondx_verified_lewis']),
    'acceptor_weak': sum(1 for r in nbo_rows if r['lewis'] == 'weak'), 'acceptor_backbone_NH': sum(1 for r in nbo_rows if r['acceptor_class'] == 'backbone N-H'),
    'verified_weak': sum(1 for r in nbo_rows if r['bondx_verified'] and r['lewis'] == 'weak'),
    'intra_min': q(min(r['E2_intra_LP'] for r in nbo_rows if r['E2_intra_LP'] > 0)), 'intra_max': q(max(r['E2_intra_LP'] for r in nbo_rows)),
    'pair_E2_range_verified_lewis': [q(min(r['E2_bondx_pair'] for r in nbo_rows if r['bondx_verified_lewis'])), q(max(r['E2_bondx_pair'] for r in nbo_rows if r['bondx_verified_lewis']))],
    'verified_lewis_list': [(r['pathway'], r['residue'], r['acceptor_atom'], r['E2_bondx_pair']) for r in nbo_rows if r['bondx_verified_lewis']],
    'verified_weak_list': [(r['pathway'], r['residue'], r['acceptor_atom'], r['E2_bondx_pair']) for r in nbo_rows if r['bondx_verified'] and r['lewis'] == 'weak']}
io = [max(x['e2'] for x in r['lpRows'] if x['frag'] == 1 and x['donor'][0] == 'O') for r in clusters if any(x['frag'] == 1 and x['donor'][0] == 'O' for x in r['lpRows'])]
orth = all(x['bondsToC'] == 2 for r in clusters for x in r['lpRows'] if x['frag'] == 1 and x['donor'][0] == 'O')
S['nbo'].update({'intra_O_clusters': len(io), 'intra_O_min': q(min(io)), 'intra_O_max': q(max(io)), 'intra_O_median': q(np.median(io)), 'intra_O_all_ortho_phenolate': orth,
                 'intra_O_ge_0_5': sum(1 for v in io if v >= 0.5)})
# ortho phenolate O term: donor in fragment 1 that is an O bonded (2 bonds away) to the C-I carbon is the dominant intramolecular term
json.dump(S, open(os.path.join(OUT, 'hiyd_summary.json'), 'w'), indent=1)
json.dump(check, open(os.path.join(OUT, 'hiyd_frame_check.json'), 'w'), indent=1)
print(json.dumps(S, indent=1))
