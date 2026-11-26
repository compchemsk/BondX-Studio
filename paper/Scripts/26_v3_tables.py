#!/usr/bin/env python3
"""26_v3_tables.py
LaTeX tables for ESI_V3 that depend on the V3 analysis (valence NBO screen with grey zone, kcal/mol units,
redesigned BXQI, competing hydrogen bonds, field tests, hIYD T2/T5). Every number comes from the result files.
Run after 11_make_tables.py (which writes the tables that did not change) and 21_v3_analysis.py.
Usage: python3 26_v3_tables.py [ROOT] [APP]
Output: tables/*.tex
"""
import os, sys, json, ast, subprocess, importlib.util
import pandas as pd, numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
APP = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, '..', 'app')
RES = os.path.join(ROOT, 'results'); V3 = os.path.join(RES, 'v3'); TAB = os.path.join(ROOT, 'tables')
spec = importlib.util.spec_from_file_location('mt', os.path.join(os.path.dirname(os.path.abspath(__file__)), '11_make_tables.py'))
mt = importlib.util.module_from_spec(spec); sys.argv = [sys.argv[0], ROOT, APP]; spec.loader.exec_module(mt)
tex, write, longtable = mt.tex, mt.write, mt.longtable
S = json.load(open(os.path.join(V3, 'summary_v3.json')))
R2 = json.load(open(os.path.join(V3, 'summary_round2.json')))
KJ = 4.184
VSHORT = {'supported': 'S', 'inconclusive': 'I', 'not supported': 'N'}
SYSN = {'1us0_idd594': '1US0', '2w0s': '2W0S', '1p5e': '1P5E', '4agl': '4AGL', '4gql': '4GQL', '4dst': '4DST'}
VARN = {'cluster_gas': 'cluster, gas', 'cluster_pcm4': r'cluster, PCM', 'shell3': r'+3~\AA{} shell', 'shell4': r'+4~\AA{} shell', 'shell5': r'+5~\AA{} shell',
        'qmmm_ff14SB': 'QM/MM ff14SB', 'qmmm_charmm36': 'QM/MM CHARMM36', 'qmmm_Z1': 'QM/MM Z1', 'qmmm_Z2': 'QM/MM Z2', 'qmmm_Z3': 'QM/MM Z3',
        'qmmm_r10': r'QM/MM 10~\AA{}', 'qmmm_r20': r'QM/MM 20~\AA{}', 'qmmm_r25': r'QM/MM 25~\AA{}', 'qmmm_full': 'QM/MM whole', 'qmmm_HIS110_HID': 'QM/MM His110 HID', 'qmmm_HIS110_HIP': 'QM/MM His110 HIP'}


def f(v, n=2):
    if v is None or (isinstance(v, float) and np.isnan(v)): return '--'
    return ('%.' + str(n) + 'f') % v


def lab(c):
    """'TBS301:A@BR11···LEU83:A@O' -> 'Br11$\\cdots$O Leu83 (A)'"""
    c = c.replace(' ··· ', '···'); d, a = c.split('···'); x = d.split('@')[1]; ch = d.split('@')[0].split(':')[1]
    r, at = a.split('@'); rn = r.split(':')[0]; res = rn[:3].title() + rn[3:]
    at = r'$\pi$' if 'π' in at else at
    return '%s$\\cdots$%s %s (%s)' % (x.title(), at, res, ch)


# ---------------------------------------------------------------- distance cut-offs
def t_cutoffs():
    B = mt.node("const B=require('%s');console.log(JSON.stringify({c:B.defaults().geom.cutoffs}))" % mt.ENG)['c']
    bondi = {'I': 1.98, 'BR': 1.85, 'CL': 1.75, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'AR': 1.77}
    rows = []
    for x in ('I', 'BR', 'CL'):
        for y in ('N', 'O', 'S', 'AR'):
            s = bondi[x] + bondi[y]; c = B[x][y]
            rows.append([x.title(), 'ring' if y == 'AR' else y, '%.2f' % c, '%.2f' % s, '%.3f' % (c / s), '%+.2f' % (c - B['I'][y]) if x != 'I' else '--'])
    body = r'''\begin{table}[!htbp]\centering\small
\caption{Default X$\cdots$Y distance cut-offs of tier T1 (\AA). The iodine values come from our earlier hIYD study.\cite{Karmakar2026} The Br and Cl values follow from them by the difference in halogen radius (I$\to$Br about $-0.13$~\AA; Br$\to$Cl about $-0.07$~\AA) and are rounded to 0.01~\AA. $\Sigma r$: sum of the Bondi radii (ring centroid: 1.77~\AA). The cut-offs lie at 1.09--1.14 times $\Sigma r$ for atom acceptors and at 1.07--1.08 times $\Sigma r$ for ring centroids. Rounding changed no T1 decision in HAL59, in the scans, or in the protein test set.}\label{tab:cutoffs}
\begin{tabular}{llrrrr}\toprule X & Y & $d_\mathrm{cut}$ & $\Sigma r$ & $d_\mathrm{cut}/\Sigma r$ & shift from I \\\midrule
''' + '\n'.join(' & '.join(r) + r'\\' for r in rows) + '\n\\bottomrule\\end{tabular}\\end{table}\n'
    write('tab_cutoffs', body)
    return [float(r[4]) for r in rows]


# ---------------------------------------------------------------- force-field reference (heavy atoms and hydrogens)
def t_ffref():
    out = subprocess.run(['node', os.path.join(APP, 'tests', 'ff_reference.test.js')], capture_output=True, text=True).stdout.split('\n')
    import re
    heavy = {}; hyd = {}
    for l in out:
        m = re.match(r'(PASS|FAIL)\s+(\S+)\s+(.+?): (\d+) residues, (\d+) heavy atoms compared, max \|Δq\| = (\S+) e(.*)', l)
        if m: heavy[(m.group(2), m.group(3))] = (m.group(4), m.group(5), float(m.group(6)), m.group(7).strip(' ()'))
        m = re.match(r'\s+(\S+)\s+(.+?): (\d+) hydrogens compared one by one, max \|Δq\| = (\S+) e(.*)', l)
        if m: hyd[(m.group(1), m.group(2))] = (m.group(3), float(m.group(4)), m.group(5).strip(' ()'))
    rows = []
    for k, (nr, na, dq, note) in heavy.items():
        h = hyd.get(k, ('--', float('nan'), ''))
        ffl = 'ff14SB vs.\\ AMBER (PDB2PQR)' if k[1].startswith('AMBER') else 'CHARMM36 vs.\\ CHARMM22 (PDB2PQR)'
        nt = 'Trp only (C36 update)' if 'CHARMM' in k[1] and dq > 0 else ''
        rows.append([k[0], ffl, nr, na, '%.2g' % dq, h[0], ('%.2g' % h[1]) if h[0] != '--' else '--', nt])
    body = r'''\begin{table}[!htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}
\caption{Force-field charges of \BX{} compared with the PDB2PQR reference files.\cite{Jurrus2018} Heavy atoms: charge of each heavy atom plus its hydrogens on both sides. Hydrogens: every hydrogen of the reference file against the template charge of the same hydrogen name, without folding. Both programs read the same AMBER parameter files, so the AMBER rows test the template matching, termini, and protonation handling of \BX{}, not the charges themselves. PDB2PQR ships CHARMM22 charges; the CHARMM36 rows only show that \BX{} differs from CHARMM22 where CHARMM36 changed the parameters (tryptophan).}\label{tab:ffref}
\resizebox{\textwidth}{!}{\begin{tabular}{llrrrrrl}\toprule PDB & Force field (BondX vs.\ reference) & Residues & Heavy atoms & max $|\Delta q|$ (e) & Hydrogens & max $|\Delta q_\mathrm{H}|$ (e) & Differences\\\midrule
''' + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\bottomrule\\end{tabular}}\\end{table}\n'
    write('tab_ffref', body)
    return sum(int(r[5]) for r in rows if r[5] != '--' and 'ff14SB' in r[1])


# ---------------------------------------------------------------- HAL59 NBO, EDA, BXQI
def t_hal59():
    q = pd.read_csv(os.path.join(RES, 'qm_hal59.csv'))
    n = pd.read_csv(os.path.join(V3, 'nbo_v3.csv')); n = n[n.set == 'A_HAL59'][['name', 'E2_valence', 'E2_RY', 'verdict', 'hbY']]
    b = pd.read_csv(os.path.join(V3, 'bxqi_hal59_table.csv'))[['name', 'BXQI', 'BXQI_old']]
    d = q.drop(columns=['BXQI', 'BXQI_grade']).merge(n, on='name').merge(b, on='name')
    rows = []
    for _, r in d.sort_values(['role', 'halogen', 'ref_kcal'], ascending=[False, True, True]).iterrows():
        role = 'control' if r.control else ('XB' if r.covered else 'XB, P/F')
        rows.append([tex(r['name']), r.halogen.title(), tex(r.acc_el), f(r.ref_kcal), f(r.G) if r.covered else '--', f(r.E2_valence), f(r.E2_RY),
                     VSHORT[r.verdict], f(r.ct_NPA * 1000, 1), f(r.eda_total_kcal), f(r.orca_cp_kcal), f(r.chi_CT), f(r.BXQI, 1), role])
    write('tab_qm_hal59', longtable('lllrrrrlrrrrrl', r'Dimer & X & Y & $E_\mathrm{ref}$ & $G$ & $\Sigma E^{(2)}$ & RY & D & $\Delta q_\mathrm{CT}$ & EDA & ORCA & $\chi_\mathrm{CT}$ & BXQI & Role', rows,
                                 r'NBO screen, ALMO-EDA, and BXQI for the 59 HAL59 dimers. $E_\mathrm{ref}$: CCSD(T)/CBS binding energy (kcal\,mol$^{-1}$; positive values mean binding). $G$: T1 score ("--": no default cut-off, tier not evaluated). $\Sigma E^{(2)}$: valence decision value (Gaussian~16, $\omega$B97X-D/def2-TZVP, kcal\,mol$^{-1}$). RY: $n\to$RY$^*$(X) part, reported beside the decision value and not counted. D: decision (S orbital-supported, $\ge 0.8$; I inconclusive, 0.3--0.8; N not supported, $<0.3$~kcal\,mol$^{-1}$). $\Delta q_\mathrm{CT}$: NPA charge moved from fragment~2 (acceptor) to fragment~1 (halogen donor), in $10^{-3}$~e. EDA: ALMO-EDA total ($\omega$B97X-V/def2-TZVPD). ORCA: counterpoise interaction energy at the same level (kcal\,mol$^{-1}$). BXQI: index from $G$, $\Sigma E^{(2)}$, and the EDA total, with tiers that cannot be evaluated left out.', 'tab:qmhal'))
    # EDA table in kcal/mol
    e = pd.read_csv(os.path.join(RES, 'qm_eda.csv')); e = e[e.set == 'A_HAL59'].merge(q[['name', 'ref_kcal', 'halogen', 'acc_el', 'orca_cp_kcal', 'role', 'control']], on='name')
    rows = []
    for _, r in e.sort_values(['role', 'halogen', 'ref_kcal'], ascending=[False, True, True]).iterrows():
        rows.append([tex(r['name']), f(r.frz / KJ), f(r.pol / KJ), f(r.ct / KJ), f(r.total / KJ), f(-r.ref_kcal), f(r.orca_cp_kcal), f(r.chi_CT)])
    write('tab_qm_eda', longtable('lrrrrrrr', r'Dimer & FRZ & POL & CT & Total & $-E_\mathrm{ref}$ & ORCA & $\chi_\mathrm{CT}$', rows,
                                r'ALMO-EDA of the 59 HAL59 dimers (Q-Chem~6.2.2, $\omega$B97X-V/def2-TZVPD, gas phase, counterpoise-corrected CT; kcal\,mol$^{-1}$). $-E_\mathrm{ref}$: CCSD(T)/CBS interaction energy. ORCA: counterpoise-corrected SCF interaction energy at the same level. With the same functional and basis, the EDA total is by construction this interaction energy; the residual difference (at most 0.05~kcal\,mol$^{-1}$) reflects RIJCOSX and integration grids.', 'tab:qmeda'))


# ---------------------------------------------------------------- protein contacts summary and per-model values
def t_protein():
    c = pd.read_csv(os.path.join(V3, 'protein_contacts_table.csv')); b = pd.read_csv(os.path.join(V3, 'bxqi_protein_table.csv'))
    c = c.merge(b[['contact', 'BXQI']], on='contact')
    rows = []
    for _, r in c.sort_values('BXQI', ascending=False).iterrows():
        v = ast.literal_eval(r.verdicts); vr = ast.literal_eval(r.verdicts_withRY_vs_05)
        fl = R2['protein_hbflag_counts'].get(r.contact.replace('···', ' ··· '), [0, r.n])
        rows.append([SYSN[r.sys], lab(r.contact), str(r.n), f(r.G), f(r.val_default), '%s--%s' % (f(r.val_min), f(r.val_max)), f(r.RY_default), '%s--%s' % (f(r.withRY_min), f(r.withRY_max)),
                     '/'.join(vr), f(r.hbY_max), f(r.hbFrag2_max), '%d/%d' % (fl[0], r.n), '/'.join(VSHORT[x] for x in v), f(r.BXQI, 1)])
    body = r'''\begin{table}[!htbp]\centering\scriptsize\setlength{\tabcolsep}{2.5pt}
\caption{Electronic tiers for the protein contacts (Gaussian~16, $\omega$B97X-D/def2-TZVP; kcal\,mol$^{-1}$). $n$: number of models. Default: PCM cluster ($\varepsilon=4$). Range: all models (gas, PCM, QM size, and QM/MM variants). RY: $n\to$RY$^*$(X) part in the default model. With RY: valence plus Rydberg sum, and the decision of a rule that counts the Rydberg terms against a single threshold of 0.5~kcal\,mol$^{-1}$ (pass/fail; the rule of BondX Studio 2.0). X$\to$Y--H: summed donation from the halogen lone pairs into the N--H and O--H antibonds of the acceptor residue; X$\to$frag.~2: the same into any N--H or O--H antibond of fragment~2. H-bond flag: models in which X$\to$frag.~2 exceeds $\Sigma E^{(2)}$ (and 0.5~kcal\,mol$^{-1}$). D: decisions over all models. BXQI: from $G$ and the default $\Sigma E^{(2)}$. All models were built with BondX Studio 2.0; the 1US0 models have the template rotamer of the Thr113 hydroxyl, whose hydrogen points at the bromine (see text).}\label{tab:protcontacts}
\resizebox{\textwidth}{!}{\begin{tabular}{llrrrlrlllrllr}\toprule PDB & Contact & $n$ & $G$ & $\Sigma E^{(2)}$ & range & RY & with RY & 0.5 rule & X$\to$Y--H & X$\to$frag.~2 & H-bond flag & D & BXQI\\\midrule
''' + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\bottomrule\\end{tabular}}\\end{table}\n'
    write('tab_protcontacts', body)
    # per model
    n = pd.read_csv(os.path.join(V3, 'nbo_v3.csv')); n = n[n.set == 'F_protein'].copy()
    n['sys'] = n.name.str.split('__').str[0]; n['var'] = n.name.str.split('__').str[-1]
    order = list(VARN)
    n['vo'] = n['var'].map(lambda v: order.index(v) if v in order else 99)
    rows = []
    for _, r in n.sort_values(['sys', 'label', 'vo']).iterrows():
        rows.append([SYSN[r.sys], lab(r.label), VARN.get(r['var'], tex(r['var'])), f(r.E2_valence), f(r.E2_RY), f(r.hbY), f(r.hbFrag2), 'yes' if str(r.lpOnY).lower() == 'true' else 'no', f(r.ct_NPA * 1000, 1), VSHORT[r.verdict]])
    write('tab_qm_prot', longtable('lllrrrrlrl', r'PDB & Contact & Model & $\Sigma E^{(2)}$ & RY & X$\to$Y--H & X$\to$frag.~2 & LP(Y) & $\Delta q_\mathrm{CT}$ & D', rows,
                                 r'NBO screen of every protein model (Gaussian~16, $\omega$B97X-D/def2-TZVP; kcal\,mol$^{-1}$). QM/MM models use ff14SB, the RC boundary, and \SI{15}{\angstrom} unless the model name says otherwise. LP(Y): NBO prints a lone pair on the acceptor atom. $\Delta q_\mathrm{CT}$: NPA charge moved to fragment~1 ($10^{-3}$~e); in a protein cluster it sums every ligand--environment contact and is not specific to the halogen bond. D: decision as in Table~\ref{tab:qmhal}.', 'tab:qmprot'))
    return len(n)


# ---------------------------------------------------------------- in-cluster pairs (hard negatives inside the models)
def t_pairs():
    p = pd.read_csv(os.path.join(V3, 'pairs_v3.csv'))
    p['sys'] = p.name.str.split('__').str[0]
    rows = []
    for _, r in p.sort_values(['pass_T1', 'E2_valence'], ascending=[True, True]).iterrows():
        why = 'T1 pass' if r.pass_T1 else ('angle' if r.pass_dist else ('distance' if r.pass_angle else 'distance, angle'))
        rows.append([SYSN[r.sys], tex(r.X.split('@')[1].title()), tex(r.Yname) + ' ' + tex(r.Yres.split(':')[0].title()) + ' (' + r.Yres.split(':')[1] + ')', f(r.d), f(r.theta1, 1), f(r.cutoff), why,
                     'poor' if str(r.lewis_poor).lower() == 'true' else 'base', f(r.E2_valence), f(r.E2_RY), f(r.hbY)])
    body = r'''\begin{table}[!htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}
\caption{Every halogen$\cdots$N, O, or S pair within 4.5~\AA{} inside the default PCM cluster models, whether or not it passes T1 (kcal\,mol$^{-1}$). The twelve pairs that pass T1 are the target contacts and their chain replicas. The six pairs that fail T1 never reach the NBO screen in the workflow, so their small values confirm the geometric exclusion rather than test the screen on contacts that pass T1; ESI Table~\ref{tab:hardneg} gives that test. Fail: criterion of T1 that is not met. Base: Lewis base or poor Lewis base. X$\to$Y--H: halogen lone pair into an N--H or O--H antibond on the partner residue.}\label{tab:pairs}
\begin{tabular}{lllrrrllrrr}\toprule PDB & X & Y & $d$ (\AA) & $\theta_1$ & $d_\mathrm{cut}$ & T1 & Base & $\Sigma E^{(2)}$ & RY & X$\to$Y--H\\\midrule
''' + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\bottomrule\\end{tabular}\\end{table}\n'
    write('tab_pairs', body)


# ---------------------------------------------------------------- MM field tests
def sgn(v):
    v = int(round(v)); return '0' if v == 0 else '%+d' % v


def t_field():
    d = pd.read_csv(os.path.join(V3, 'field_v3.csv'))
    rows = []
    for s in SYSN:
        g = d[d.system == s]
        def E(ff='ff14SB', bd='RC', rad='15', eh=0):
            q = g[(g.ff == ff) & (g.boundary == bd) & (g.radius.astype(str) == rad) & (g.explicitH == eh)]
            return q.iloc[0] if len(q) else None
        a, b, c, w, ch = E(eh=0), E(eh=6), E(bd='RCD', eh=6), E(rad='full', eh=6), E(ff='charmm36', eh=6)
        rows.append([SYSN[s], f(a.E_AX), f(b.E_AX), f(c.E_AX), f(ch.E_AX), f(E(rad='25', eh=6).E_AX), f(w.E_AX), '%d' % b.nExplicitH, '%d/%d' % (b.nZero, ch.nZero), '%d/%d' % (b.nCharges, ch.nCharges), sgn(b.net), sgn(w.net)])
    body = r'''\begin{table}[!htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}
\caption{Field component $E_{AX}$ (V\,\AA$^{-1}$) along the A$\to$X bond at the halogen from the MM point charges of the QM/MM models (no QM job), with the QM regions of the computed models. The hydrogen-bond shell of version 3.0 moves Lys105 of 2W0S into the QM region, which changes $E_{AX}$ from the MM charges of the whole structure from $-0.44$ to $+0.22$~V\,\AA$^{-1}$; the field values therefore refer to a given QM region. Positive values point from A to X and should deepen the $\sigma$-hole. Folded: hydrogen charges moved onto their heavy atoms everywhere (the setting of the QM/MM models of this work, built with BondX Studio 2.0). Explicit H: MM residues within 6~\AA{} of X carry explicit hydrogens with their template charges (default of version 3.0). $n_\mathrm{H}$: explicit MM hydrogens. Zero: MM atoms with zero charge after folding (ff14SB/CHARMM36); CHARMM36 gives zero net charge to many aliphatic CH$_n$ groups, so after folding these atoms carry no charge and are not written, which explains the smaller CHARMM36 point-charge counts at the same radius. $N_q$: written point charges (ff14SB/CHARMM36, \SI{15}{\angstrom}). $Q$: MM net charge (e) at \SI{15}{\angstrom} and for the whole structure.}\label{tab:field}
\resizebox{\textwidth}{!}{\begin{tabular}{lrrrrrrrrrrr}\toprule PDB & folded & explicit H & RCD & CHARMM36 & 25~\AA & whole & $n_\mathrm{H}$ & Zero & $N_q$ & $Q_{15}$ & $Q_\mathrm{whole}$\\\midrule
''' + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\bottomrule\\end{tabular}}\\end{table}\n'
    write('tab_field', body)


def t_m12():
    m = S['m12']; rows = []
    W = {w[0]: w for w in m['within']}
    for c, e, dv, dp in m['rows_default']:
        w = W.get(c, [c, float('nan'), 0])
        rows.append([lab(c), f(e), '%+.2f' % dv, '%+.2f' % dp, f(w[1]) if w[1] == w[1] else '--', str(w[2])])
    body = r'''\begin{table}[!htbp]\centering\small
\caption{Does the MM field predict the change of $\Sigma E^{(2)}$ on embedding? $E_{AX}$: default QM/MM model (V\,\AA$^{-1}$). $\Delta\Sigma E^{(2)}$: QM/MM minus gas-phase cluster, and QM/MM minus PCM cluster (kcal\,mol$^{-1}$). $\rho$: Spearman correlation between $E_{AX}$ and $\Sigma E^{(2)}$ over the QM/MM variants of one contact ($n$ models). Across contacts, $\rho = %.2f$ ($n = %d$); over all variants pooled, $\rho = %.2f$ ($n = %d$).}\label{tab:m12}
\begin{tabular}{lrrrrr}\toprule Contact & $E_{AX}$ & vs.\ gas & vs.\ PCM & $\rho$ within & $n$\\\midrule
''' % (m['rho_across_default'], m['n_contacts'], m['rho_all_variants'], m['n_all']) + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\bottomrule\\end{tabular}\\end{table}\n'
    write('tab_m12', body)


def t_levels():
    L = S['levels']['rel']; rows = []
    for k in ['HF__def2-TZVP', 'B3LYP-D3BJ__def2-TZVP', 'PBE0-D3BJ__def2-TZVP', 'M06-2X__def2-TZVP', 'wB97X-D__def2-SVP', 'wB97X-D__def2-QZVP']:
        m, b = k.split('__'); m = m.replace('wB97X-D', r'$\omega$B97X-D').replace('D3BJ', 'D3(BJ)')
        rows.append([m, b, f(L[k]['min']), f(L[k]['max'])])
    body = r'''\begin{table}[!htbp]\centering\small
\caption{Valence $\Sigma E^{(2)}$ of the 12 scan dimers at $d_\mathrm{eq}$ relative to $\omega$B97X-D/def2-TZVP (Set D, Gaussian~16). No decision changes; the smallest value at any level is %.2f~kcal\,mol$^{-1}$; the rank correlation with the default level is at least %.2f.}\label{tab:levels}
\begin{tabular}{llrr}\toprule Method & Basis & min ratio & max ratio\\\midrule
''' % (S['levels']['min_val'], S['levels']['rho_min']) + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\bottomrule\\end{tabular}\\end{table}\n'
    write('tab_levels', body)


def t_protein_eda():
    rows = []
    for s, lig, frz, pol, ct, tot, chi, g16, orca in S['protein_eda']:
        rows.append([SYSN[s], tex(lig[:-1]) + ' (' + lig[-1] + ')', f(frz), f(pol), f(ct), f(tot), f(chi), f(orca), f(g16)])
    body = r'''\begin{table}[!htbp]\centering\small
\caption{ALMO-EDA of the six finished protein-cluster jobs (Q-Chem~6.2.2, $\omega$B97X-V/def2-TZVPD, gas phase; kcal\,mol$^{-1}$) with the counterpoise interaction energies from ORCA (same level) and Gaussian~16 ($\omega$B97X-D/def2-TZVP, no diffuse functions). The terms cover every contact between the ligand and the acceptor fragment. The large Gaussian--ORCA gap for 2W0S (ligand charge $-2$) shows the need for diffuse functions on anions. Three EDA jobs (2W0S chain A, 1P5E chain A, and 4GQL) stopped during the fragment calculations and gave no result; the 4GQL ORCA job was not run.}\label{tab:proteda}
\begin{tabular}{llrrrrrrr}\toprule PDB & Ligand & FRZ & POL & CT & Total & $\chi_\mathrm{CT}$ & ORCA & G16\\\midrule
''' + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\bottomrule\\end{tabular}\\end{table}\n'
    write('tab_proteda', body)


# ---------------------------------------------------------------- hIYD NBO and T2/T5
def t_hiyd():
    h = pd.read_csv(os.path.join(V3, 'hiyd_decisions_updated.csv'))
    rows = []
    for _, r in h.sort_values(['pathway', 'xi']).iterrows():
        bb = r.source.startswith('acetyl')
        rows.append([r.pathway, tex(r.residue), f(r.xi, 1), tex(r.acceptor_atom), 'poor' if r.lewis == 'weak' else 'base', f(r.d_A), f(r.theta1, 1), f(r.E_pub), f(r.E2_intra_phenolateO),
                     f(r.E2_bondx_pair), VSHORT[r.verdict_v3], f(r.E2_used) if bb else '', f(r.hb_used), VSHORT[r.verdict_used] + (r'$^\dagger$' if r.hbFlag else '')])
    write('tab_hiyd_nbo', longtable('llrllrrrrrlrrl', r'Path & Residue & $\xi$ & Y & Base & $d$ & $\theta_1$ & Pub. & Intra & $\Sigma E^{(2)}$ & D & cap & X$\to$Y--H & Final', rows,
                                  r'NBO screen of the 56 published hIYD poses (B3LYP-D3/6-311++G(d,p) with LANL2DZdp on I, CPCM $\varepsilon=4$; kcal\,mol$^{-1}$).\cite{Karmakar2026} $\xi$: window (\AA). $d$ (\AA) and $\theta_1$ ($^\circ$): T1 geometry of the published pose. Base: Lewis class of the acceptor in the protein (poor: delocalized or protonated lone pair). Pub.: published value. Intra: lone pair of the phenolate O of I-Tyr into $\sigma^*$(C--I) of the same ring. $\Sigma E^{(2)}$ and D: valence decision value and decision for the published model; all 18 backbone N acceptors of the published models are amines, because a hydrogen replaced the preceding carbonyl carbon. cap: valence value after this hydrogen was rebuilt as an acetyl group (set K; same pose and level). X$\to$Y--H: halogen lone pairs into N--H or O--H antibonds (for backbone N, of the acetyl-capped model). Final: decision with the acetyl-capped models for backbone N; $^\dagger$: the competing hydrogen bond to I exceeds $\Sigma E^{(2)}$. One set-K job (A130$^*$, pathway C) stopped while writing the NBO file after the complete NBO analysis had been printed; its values are used. These outputs use a different level of theory from the one the grey zone was derived for.', 'tab:hiydnbo', size=r'\scriptsize'))
    t = pd.read_csv(os.path.join(V3, 'hiyd_t5_table.csv'))
    import re as _re
    # the stripped trajectories number residues one higher than the published cluster files (see 21_v3_analysis.py)
    amap = {(_re.match(r'([A-Z]+)(\d+)', a).group(1) + str(int(_re.match(r'([A-Z]+)(\d+)', a).group(2)) + 1)): r for a, r in zip(h.residue_amber, h.residue)}
    upd = {(x.pathway, x.residue, x.acceptor_atom, round(x.xi, 1)): x.E2_used for _, x in h.iterrows()}
    rows = []
    for _, r in t.iterrows():
        p, res, at = r.key.split('|')
        same = (r.xi_pub == r.xi_pub) and abs(r.xi_pub - r.xi) < 1e-6 and (r.E2_bondx_pair == r.E2_bondx_pair)
        rows.append([p, tex(amap.get(res, res.title()[:3] + res[3:])), at, f(r.xi, 1), f(r.occupancy * 100, 1), '%.0f' % r.tau_ps, f(r.D), f(r.D / r.occupancy), f(r.G_best),
                     f(upd.get((p, amap.get(res, ''), at, round(r.xi, 1)), r.E2_bondx_pair)) if same else '--'])
    body = r"""\begin{table}[!htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}
\caption{Tier T2 on the hIYD release windows for the contacts with the highest occupancy. Occupancy: share of the frames of one umbrella window (unweighted, \%). $\tau$: mean residence time (ps; 20~ps per MD frame times the stride of the stripped file). $D$: persistence score ($\tau_0 = 100$~ps); at the frame spacing of these windows (160--240~ps), $D/f \ge 0.95$, so $D$ is practically the occupancy. $G$: best frame. $\Sigma E^{(2)}$: valence value of the NBO output for the same residue and acceptor atom in the same window (published model; acetyl-capped model for backbone N; ``--'': no output from this window).}\label{tab:hiydt5}
\begin{tabular}{lllrrrrrrr}\toprule Path & Residue & Y & $\xi$ & Occ. & $\tau$ & $D$ & $D/f$ & $G$ & $\Sigma E^{(2)}$\\\midrule
""" + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\bottomrule\\end{tabular}\\end{table}\n'
    write('tab_hiydt5', body)


# ---------------------------------------------------------------- round-2 tables
def t_classratio():
    c = pd.read_csv(os.path.join(V3, 'class_ratio.csv'), index_col=0)
    nm = {'P': 'P', 'N': 'N', 'S': 'S', 'O': 'O', 'F': 'F', 'ring': r'ring ($\pi$)'}
    rows = [[nm.get(k, k), '%d' % v['count'], f(v['median']), '%s--%s' % (f(v['min']), f(v['max'])), f(0.8 / v['median'])] for k, v in c.iterrows()]
    body = r"""\begin{table}[!htbp]\centering\small
\caption{Valence $\Sigma E^{(2)}$ per unit CCSD(T)/CBS binding energy for each acceptor class of the 56 HAL59 halogen bonds. BE$_{0.8}$: binding energy (kcal\,mol$^{-1}$) that corresponds to the upper grey-zone limit of 0.8~kcal\,mol$^{-1}$ at the median ratio. Ring acceptors are under-scored about fourfold relative to N acceptors; with only four ring dimers, no class-specific zone is fitted.}\label{tab:classratio}
\begin{tabular}{lrrlr}\toprule Acceptor & $n$ & median & range & BE$_{0.8}$\\\midrule
""" + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\bottomrule\\end{tabular}\\end{table}\n'
    write('tab_classratio', body)


def t_hardneg():
    ex = pd.read_csv(os.path.join(V3, 'extra_nbo.csv')); ex['flags'] = ex['flags'].fillna('')
    G = ex[ex.set == 'G_hard_negatives'].copy()
    G['kind'] = np.where(G.role.str.startswith('control'), 'base', 'poor')
    G['res'] = G.role.str.extract(r'\(([^)]+)\)')[0]; G['what'] = G.role.str.replace(r'^(negative|control): ', '', regex=True).str.replace(r' \([^)]+\)', '', regex=True)
    G['path'] = G.name.str.extract(r'hiyd_C(\d)_')[0].map({'1': 'A', '2': 'B', '3': 'C'}); G['xi'] = G.name.str.extract(r'hiyd_C\d_([\d.]+)_')[0].astype(float)
    rows = []
    for _, r in G.sort_values(['kind', 'E2_valence'], ascending=[False, True]).iterrows():
        rows.append([tex(r.res), tex(r.what), r.kind, r.path, f(r.xi, 1), '%d' % r.frame, f(r.d), f(r.theta1, 1), f(r.G), f(r.E2_valence), f(r.E2_RY), f(max(r.hbY, r.hbFrag2)), 'yes' if str(r.lpOnY).lower() == 'true' else 'no', VSHORT[r.verdict] + (r'$^\dagger$' if 'hbond' in r['flags'] else '')])
    J = ex[ex.set == 'J_4dst_amide']
    ENVN = {'cluster_gas': 'gas', 'cluster_pcm4': 'PCM', 'qmmm_ff14SB': 'QM/MM'}
    jr = [['Gly75 (4DST, %s)' % ENVN[r['name'].split('__')[2]], 'backbone amide N', 'poor', '--', '--', '--', '3.52', '140.5', '0.12', f(r.E2_valence), f(r.E2_RY), f(max(r.hbY, r.hbFrag2)), 'yes' if str(r.lpOnY).lower() == 'true' else 'no', VSHORT[r.verdict]] for _, r in J.iterrows()]
    envs = [n.split('__')[2] for n in J.name]
    body = r"""\begin{table}[!htbp]\centering\scriptsize\setlength{\tabcolsep}{2.6pt}
\caption{Contacts that pass T1 but involve poor Lewis bases (hard negatives), with Lewis-base controls from the same simulations (kcal\,mol$^{-1}$). hIYD rows: one frame of a release window each, cluster models built with the defaults of BondX Studio 3.0 (intact amides, oriented hydroxyl hydrogens), I-Tyr phenolate charge $-1$, CPCM $\varepsilon=4$, $\omega$B97X-D/def2-TZVP (Gaussian~16, NBO~3.1); every output matched its manifest. 4DST rows: the Cl$\cdots$N Gly75 contact with the amide of Gly75 intact (""" + ', '.join(envs).replace('cluster_gas', 'gas cluster').replace('cluster_pcm4', 'PCM cluster').replace('qmmm_ff14SB', 'QM/MM') + r"""). Base: Lewis base or poor Lewis base. X$\to$H: summed donation from the iodine lone pairs into the N--H and O--H antibonds of fragment~2. D: decision; $^\dagger$: the competing hydrogen bond exceeds $\Sigma E^{(2)}$.}\label{tab:hardneg}
\resizebox{\textwidth}{!}{\begin{tabular}{lllcrrrrrrrrll}\toprule Residue & Acceptor & Base & Path & $\xi$ & Frame & $d$ & $\theta_1$ & $G$ & $\Sigma E^{(2)}$ & RY & X$\to$H & LP(Y) & D\\\midrule
""" + '\n'.join(' & '.join(x) + r'\\' for x in rows) + '\n\\midrule\n' + '\n'.join(' & '.join(x) + r'\\' for x in jr) + '\n\\bottomrule\\end{tabular}}\\end{table}\n'
    write('tab_hardneg', body)


if __name__ == '__main__':
    print('cut-off ratios', t_cutoffs())
    print('hydrogens compared (AMBER)', t_ffref())
    t_hal59(); print('protein NBO rows', t_protein()); t_pairs(); t_field(); t_m12(); t_levels(); t_protein_eda(); t_hiyd(); t_classratio(); t_hardneg()
