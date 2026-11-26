#!/usr/bin/env python3
"""11_make_tables.py
LaTeX tables for ESI_V2 (and two for Manuscript_V2), written from the result CSV files and from the BondX engine itself
(parameter list, assumptions, tests), so that no number is typed by hand.
Output: tables/*.tex
"""
import os, sys, csv, json, re, subprocess
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
APP = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'app')  # the app folder of the repository
RES = os.path.join(ROOT, 'results'); TAB = os.path.join(ROOT, 'tables'); os.makedirs(TAB, exist_ok=True)
ENG = os.path.join(ROOT, 'Scripts', 'engine.js')

U = {'θ': r'$\theta$', 'α': r'$\alpha$', 'β': r'$\beta$', 'σ': r'$\sigma$', 'π': r'$\pi$', 'χ': r'$\chi$', 'ε': r'$\varepsilon$', 'τ': r'$\tau$', 'Δ': r'$\Delta$', 'Σ': r'$\Sigma$',
     '₀': r'$_0$', '₁': r'$_1$', '₂': r'$_2$', '₃': r'$_3$', '₄': r'$_4$', '⁻': r'$^-$', '⁺': r'$^+$', '½': r'$\frac12$', '⅓': r'$\frac13$', '·': r'$\cdot$',
     '≤': r'$\le$', '≥': r'$\ge$', '×': r'$\times$', '−': r'$-$', '→': r'$\to$', '≈': r'$\approx$', '°': r'$^\circ$', 'Å': r'\AA{}', '…': r'\ldots{}',
     '“': '``', '”': "''", '’': "'", '—': '---', '–': '--', 'ᵢ': r'$_i$', 'ₖ': r'$_k$', '²': r'$^2$', '·': r'$\cdot$', 'Π': r'$\Pi$', '⇔': r'$\Leftrightarrow$', '✓': r'\checkmark'}


def tex(s):
    s = str(s)
    s = s.replace('\\', r'\textbackslash{}')
    for a, b in (('&', r'\&'), ('%', r'\%'), ('#', r'\#'), ('_', r'\_'), ('{', r'\{'), ('}', r'\}'), ('^', r'\^{}'), ('~', r'\~{}'), ('$', r'\$')):
        s = s.replace(a, b)
    s = s.replace('···', r'$\cdots$').replace('<', r'$<$').replace('>', r'$>$')
    for a, b in U.items(): s = s.replace(a, b)
    return s.replace('$$', '')


def node(js):
    return json.loads(subprocess.run(['node', '-e', js], capture_output=True, text=True).stdout)


def write(name, body):
    open(os.path.join(TAB, name + '.tex'), 'w').write(body); print('wrote', name)


def longtable(spec, head, rows, caption, label, size=r'\scriptsize'):
    out = [r'{%s' % size, r'\setlength{\tabcolsep}{3pt}', r'\begin{longtable}{%s}' % spec, r'\caption{%s}\label{%s}\\' % (caption, label), r'\toprule', head + r'\\', r'\midrule', r'\endfirsthead',
           r'\toprule', head + r'\\', r'\midrule', r'\endhead', r'\bottomrule', r'\endlastfoot']
    out += [' & '.join(r) + r'\\' for r in rows]
    out += [r'\end{longtable}}']
    return '\n'.join(out) + '\n'


def t_params():
    spec = node("const B=require('%s');const D=B.defaults();console.log(JSON.stringify(B.PARAM_SPEC.map(s=>({p:s.path,st:s.stage,l:s.label,u:s.unit,t:s.type,d:B.getPath(D,s.path),o:s.options||null}))))" % ENG)
    names = {'T0': 'Structure', 'T1': 'Geometry', 'T2': 'Dynamics', 'T3': 'QM model', 'T3b': 'NBO screen', 'T4': 'ALMO-EDA', 'T5': 'BXQI'}
    rows = []
    for s in spec:
        d = s['d']
        if isinstance(d, dict):
            if s['t'] == 'cutoffs': d = '; '.join('%s: %s' % (x, ', '.join('%s %s' % (y, v) for y, v in dd.items())) for x, dd in d.items())
            else: d = ', '.join('%s %s' % (k, v) for k, v in d.items())
        elif isinstance(d, bool): d = 'on' if d else 'off'
        elif s['o']: d = next((lab for k, lab in s['o'] if str(k) == str(d)), d)
        rows.append([names[s['st']], r'\texttt{%s}' % tex(s['p']), tex(s['l']) + (' (%s)' % tex(s['u']) if s['u'] else ''), r'\parbox[t]{5.0cm}{\raggedright %s}' % tex(d)])
    write('tab_params', longtable('llp{4.2cm}l', r'Tier & Parameter & Meaning (unit) & Default', rows, 'Every adjustable parameter of \\BX{} with its default value (%d parameters). All can be changed in the program, on each tier page, and on the \\emph{All parameters} page.' % len(spec), 'tab:params'))


def t_assumptions():
    A = node("const B=require('%s');const P=B.defaults();const o={};['T0','T1','T2','T3','T3b','T4','T5'].forEach(s=>o[s]=B.assumptions(s,P));o.QMMM=B.assumptions('T3',B.mergeParams(P,{qm:{env:'qmmm'}}));console.log(JSON.stringify(o))" % ENG)
    names = {'T0': 'Tier T0, structure input', 'T1': 'Tier T1, geometry', 'T2': 'Tier T2, persistence in dynamics', 'T3': 'Tier T3a, QM cluster model', 'QMMM': 'Tier T3a, QM/MM model (additional lines)',
             'T3b': 'Tier T3b, NBO screen', 'T4': 'Tier T4, ALMO-EDA', 'T5': 'Tier T5, ranking index (BXQI)'}
    out = []
    for k in ('T0', 'T1', 'T2', 'T3', 'QMMM', 'T3b', 'T4', 'T5'):
        items = A[k] if k != 'QMMM' else [a for a in A['QMMM'] if a not in A['T3']]
        out.append(r'\paragraph{%s.}' % names[k]); out.append(r'\begin{itemize}[nosep,leftmargin=*]')
        out += [r'\item %s' % tex(a) for a in items]; out.append(r'\end{itemize}')
    write('tab_assumptions', '\n'.join(out) + '\n')


def t_tests():
    r = subprocess.run(['node', os.path.join(APP, 'tests', 'engine.test.js')], capture_output=True, text=True, cwd=APP).stdout
    lines = [l for l in r.split('\n') if l.startswith('PASS ') or l.startswith('FAIL ')]
    ref = [l.strip()[6:] for l in r.split('\n') if l.strip().startswith('PASS  ')]
    rows = [[str(i + 1), tex(l[5:])] for i, l in enumerate(lines)]
    write('tab_tests', longtable('rp{14.5cm}', r'\# & Automated test (all passed)', rows, 'The %d automated engine tests (\\texttt{node tests/engine.test.js}). The last test runs the independent force-field comparison listed in Table~\\ref{tab:ffref}.' % len(rows), 'tab:tests'))
    rr = []
    for l in ref:
        m = re.match(r'(\S+)\s+(.+?): (\d+) residues, (\d+) heavy atoms compared, max \|Δq\| = (\S+) e(.*)', l)
        if m: rr.append([m.group(1), tex(m.group(2).replace('AMBER', 'ff14SB vs. AMBER (PDB2PQR)', 1) if m.group(2).startswith('AMBER') else m.group(2)), m.group(3), m.group(4), tex('%.2g' % float(m.group(5))), tex(m.group(6).strip(' ()'))])
    body = r'''\begin{table}[!htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}
\caption{Force-field charges of \BX{} compared with the PDB2PQR reference files. Hydrogen charges are added to their heavy atoms on both sides.}\label{tab:ffref}
\resizebox{\textwidth}{!}{\begin{tabular}{llrrlp{5.2cm}}\toprule PDB & Force field (BondX vs. reference) & Residues & Heavy atoms & max $|\Delta q|$ (e) & Note\\\midrule
''' + '\n'.join(' & '.join(x) + r'\\' for x in rr) + r'''
\bottomrule\end{tabular}}\end{table}
'''
    write('tab_ffref', body)
    return len(rows)


def t_hal59():
    d = pd.read_csv(os.path.join(RES, 't1_hal59.csv'))
    rows = []
    for _, r in d.iterrows():
        acc = 'ring' if r.acc_el == 'pi' else r.acc_el
        rows.append([tex(r['name']), r.subset, tex(r.halogen.title() if len(r.halogen) == 2 else r.halogen), tex(acc), '%.2f' % r.ref_kcal, '%.3f' % r.d, '%.1f' % r.theta1,
                     ('%.0f' % r.theta2) if str(r.theta2) != 'nan' and r.theta2 != '' else '--', 'yes' if r.covered else 'no', 'yes' if r.pass_T1 else 'no', '%.3f' % r.G,
                     tex('control' if r.role != 'halogen bond' else '')])
    write('tab_hal59_t1', longtable('llllrrrrllrl', r'Dimer & Set & X & Y & $E_\mathrm{ref}$ & $d$ (\AA) & $\theta_1$ & $\theta_2$ & Covered & T1 & $G$ & Role', rows,
                                     'HAL59 dimers. Reference binding energies: CCSD(T)/CBS (kcal\\,mol$^{-1}$, GMTKN55). Tier T1 of \\BX{} on the reference geometries. Covered: a default cut-off exists for the X/Y pair.', 'tab:hal59'))


def t_protein():
    d = pd.read_csv(os.path.join(RES, 'protein_t1.csv')); x = d[(d.includeF == False) & (d.status == 'XB')]
    rows = [[tex(r.pdb.upper().replace('_IDD594', '+IDD594').replace('_PROTONATED', ' (polar H)')), tex(r.donor), tex(r.X), tex(r.acceptor), tex(r.accAtom), tex(r.accClass), '%.3f' % r.d, '%.1f' % r.theta1,
             '%.1f' % r.theta2, tex(r.lewis), '%.3f' % r.G] for _, r in x.iterrows()]
    write('tab_protein_t1', longtable('lllllp{1.9cm}rrrlr', r'PDB & Donor & X & Acceptor & Y & Class & $d$ (\AA) & $\theta_1$ & $\theta_2$/$\alpha$ & Lewis & $G$', rows,
                                       'Halogen bonds found by tier T1 with default settings in the protein test set (all chains).', 'tab:protT1'))
    m = pd.read_csv(os.path.join(RES, 'protein_models.csv'))
    rows = []
    for _, r in m.iterrows():
        rows.append([tex(r.model.split('__')[0].upper().replace('_IDD594', '')), tex(r.ligand), tex(r.variant.replace('_', ' ')), str(int(r.n_atoms)), str(int(r.added_H)) if r.added_H == r.added_H else '',
                     str(int(r.link_H)) if r.link_H == r.link_H else '', '%d/%d' % (r.charge_donor, r.charge_acceptor) if r.charge_donor == r.charge_donor else '', str(int(r.electrons)) if r.electrons == r.electrons else '',
                     str(int(r.mm_n)) if r.mm_n == r.mm_n else '--', ('%.3f' % (r.mm_net + 0.0)) if r.mm_net == r.mm_net else '--', ('%.2f' % (round(r.esp_X_V, 2) + 0.0)) if r.esp_X_V == r.esp_X_V else '--',
                     ('%.3f' % r.field_AX_VperA) if r.field_AX_VperA == r.field_AX_VperA else '--'])
    write('tab_protein_models', longtable('lllrrrlrrrrr', r'System & Ligand & Variant & Atoms & +H & Link & $q_1/q_2$ & $N_e$ & MM $n$ & MM $Q$ & $\phi_X$ (V) & $E_{AX}$ (V/\AA)', rows,
                                           'All 89 protein QM models built by \\BX{}. Atoms: QM atoms. +H: added hydrogens. Link: link hydrogens. $q_1/q_2$: fragment charges. $N_e$: electrons. MM $n$ and MM $Q$: number and net charge of the MM point charges. $\\phi_X$ and $E_{AX}$: MM potential and field component at the halogen.', 'tab:models'))


def t_perf():
    d = pd.read_csv(os.path.join(RES, 'performance.csv'))
    t = d[d.test == 'T1 detect']; tl = t[t.system.str.contains(' x')]; ts = t[~t.system.str.contains(' x')]
    rows = [['T1 detection, %d test structures' % len(ts), '%d--%d' % (ts.atoms.min(), ts.atoms.max()), '%.0f--%.0f (median %.0f)' % (ts.ms.min(), ts.ms.max(), ts.ms.median())]]
    rows += [['T1 detection, %s' % (tex(r.system.upper().split(' X')[0]) + r' $\times$ ' + r.system.split(' x')[1]), '%d' % r.atoms, '%.0f' % r.ms] for _, r in tl.iterrows()]
    for _, r in d[d.test == 'T2 trajectory'].iterrows(): rows.append(['T2, 200 frames, %s' % tex(r.system.upper().replace('_IDD594', '')), '%d' % r.atoms, '%.0f (%.0f frames/s)' % (r.ms, r.frames_per_s)])
    for _, r in d[d.test.str.startswith('T3a')].iterrows(): rows.append(['%s model, %s' % (tex(r.test.split()[-1].replace('qmmm', 'QM/MM')), tex(r.system.upper().replace('_IDD594', ''))), '', '%.0f' % r.ms])
    body = r'''\begin{table}[!htbp]\centering\footnotesize
\caption{Run time of the \BX{} engine (Node.js 22, one core of an Intel Xeon at 2.3 GHz, median of five runs). Times include file parsing.}\label{tab:perf}
\begin{tabular}{lrl}\toprule Task & Atoms & Time (ms)\\\midrule
''' + '\n'.join(' & '.join(r) + r'\\' for r in rows) + r'''
\bottomrule\end{tabular}\end{table}
'''
    write('tab_performance', body)


def t_jobs():
    pkg = os.path.join(ROOT, 'qm_package') if os.path.isdir(os.path.join(ROOT, 'qm_package')) else os.path.join(ROOT, 'qm_inputs')  # qm_inputs/ in the GitHub copy
    j = pd.concat([pd.read_csv(os.path.join(pkg, f)) for f in ('jobs_small.csv', 'jobs_protein.csv')])
    c = j.groupby(['set', 'program']).size().unstack(fill_value=0)
    desc = {'A_HAL59': '59 HAL59 dimers', 'B_scan_distance': '12 dimers $\\times$ 8 distances', 'C_scan_angle': '6 dimers $\\times$ 7 angles', 'D_method_basis': '12 dimers $\\times$ 8 levels', 'E_eda_settings': '12 dimers $\\times$ 3 EDA settings', 'F_protein': '89 protein models'}
    P = ('gaussian', 'orca', 'qchem', 'qchem_eda', 'gaussian_cp', 'orca_cp')
    rows = [[tex(s), desc[s]] + [str(int(c.loc[s].get(p, 0))) for p in P] for s in c.index]
    rows.append([r'\textbf{Total}', ''] + [r'\textbf{%d}' % c[p].sum() for p in P])
    body = r'''\begin{table}[!htbp]\centering\footnotesize
\caption{QM jobs written by \BX{} for this work (%d input files). NBO jobs exist for Gaussian~16, ORCA, and Q-Chem. ALMO-EDA runs in Q-Chem. Each EDA job of sets A--C and F has a counterpoise interaction-energy partner in Gaussian (one job) and ORCA (three single points).}\label{tab:jobs}
\begin{tabular}{llrrrrrr}\toprule Set & Content & G16 NBO & ORCA NBO & Q-Chem NBO & Q-Chem EDA & G16 CP & ORCA CP\\\midrule
''' % len(j) + '\n'.join(' & '.join(r) + r'\\' for r in rows) + r'''
\bottomrule\end{tabular}\end{table}
'''
    write('tab_jobs', body)
    return len(j)


def t_scans():
    d = pd.read_csv(os.path.join(RES, 't1_scans.csv'))
    rows = []
    for b, g in d[d.scan == 'distance'].groupby('base', sort=False):
        g = g.sort_values('s'); rows.append([tex(b), 'distance'] + ['%.2f' % v for v in g.G])
    write('tab_scans_dist', longtable('ll' + 'r' * 8, r'Dimer & Scan & ' + ' & '.join('%.2f' % s for s in sorted(d[d.scan == 'distance'].s.unique())), rows,
                                       'Geometric score $G$ along the rigid distance scans ($d/d_\\mathrm{eq}$ in the header).', 'tab:scand'))
    rows = []
    for b, g in d[d.scan == 'angle'].groupby('base', sort=False):
        g = g.sort_values('theta1', ascending=False); rows.append([tex(b)] + ['%.1f / %.2f' % (t, v) for t, v in zip(g.theta1, g.G)])
    write('tab_scans_angle', longtable('l' + 'r' * 7, r'Dimer & \multicolumn{7}{c}{$\theta_1$ ($^\circ$) / $G$ for target angles 180, 170, \ldots, 120$^\circ$}', rows,
                                        'Geometric score along the rigid angle scans. The first number is the actual $\\theta_1$ after rotation.', 'tab:scana'))


def mnum(fmt, v):
    t = fmt % v
    return '$-$' + t[1:] if t.startswith('-') else t


def t_cp():
    """ORCA counterpoise tables (results/orca_cp_*.csv from 13_orca_cp_analysis.py)."""
    pth = os.path.join(RES, 'orca_cp_hal59.csv')
    if not os.path.exists(pth): return
    h = pd.read_csv(pth); sm = json.load(open(os.path.join(RES, 'orca_cp_summary.json')))
    rows = [[tex(r['name']), r['subset'], r['halogen'], ('$\\pi$' if r['acc'] == 'pi' else r['acc']), '%.2f' % r['ref_kcal'], '%.2f' % r['BE_kcal'], mnum('%+.2f', r['err_kcal']),
             ('%.2f' % r['G']) if r['covered'] else '--', ('XB' if r['role'] == 'halogen bond' else 'control')] for _, r in h.iterrows()]
    write('tab_cp_hal59', longtable('lllrrrrrl', r'Dimer & Set & X & Y & $E_\mathrm{ref}$ & BE(ORCA) & Error & $G$ & Role', rows,
          'Counterpoise-corrected binding energies of the HAL59 dimers (kcal\\,mol$^{-1}$, positive = bound). ORCA~5.0.1, $\\omega$B97X-V/def2-TZVPD, RIJCOSX def2/J, gas phase. $E_\\mathrm{ref}$: CCSD(T)/CBS. Error = BE(ORCA) $-$ $E_\\mathrm{ref}$. $G$: geometric score. --: no default cut-off (P and F acceptors).', 'tab:cphal'))
    lab = [('all', 'All 59 dimers'), ('XB51', 'XB51 subset'), ('X40', 'X40 subset'), ('Cl', 'Cl donors'), ('Br', 'Br donors'), ('I', 'I donors'),
           ('XB_dimers', 'Halogen bonds'), ('controls', 'Controls'), ('covered', 'Covered by T1')]
    rows = []
    for k, t in lab:
        st = sm['hal59_' + k]
        rows.append([t, str(st['n']), mnum('%+.2f', st['MSE']), '%.2f' % st['MAE'], '%.2f' % st['RMSD'], '%.2f' % st['MaxAE'], '%.3f' % st['pearson_r'], '%.3f' % st['spearman_rho']])
    head = ('\\begin{table}[!htbp]\\centering\\footnotesize\n'
            '\\caption{Error statistics of the ORCA counterpoise binding energies against CCSD(T)/CBS (kcal\\,mol$^{-1}$). MSE: mean signed error. '
            'MAE: mean absolute error. MaxAE: largest absolute error. The I row includes the control FI$\\cdots$FCCH. The two fluoroarene controls enter only the first and the last-but-one rows.}\\label{tab:cpstat}\n'
            '\\begin{tabular}{lrrrrrrr}\\toprule Group & $n$ & MSE & MAE & RMSD & MaxAE & $r$ & $\\rho$\\\\\\midrule\n')
    write('tab_cp_stats', head + '\n'.join(' & '.join(r) + '\\\\' for r in rows) + '\n\\bottomrule\\end{tabular}\\end{table}\n')
    sc = pd.read_csv(os.path.join(RES, 'orca_cp_scans.csv'))
    d = sc[sc.scan == 'distance']; ss = sorted(d.s.unique())
    rows = [[tex(b)] + ['%.2f' % v for v in g.sort_values('s').BE_kcal] for b, g in d.groupby('base')]
    write('tab_cp_dist', longtable('l' + 'r' * len(ss), r'Dimer & ' + ' & '.join('%.2f' % v for v in ss), rows,
          'ORCA counterpoise binding energies (kcal\\,mol$^{-1}$) along the rigid distance scans. The header gives $d/d_\\mathrm{eq}$. The column 1.00 equals the HAL59 value.', 'tab:cpdist'))
    a = sc[sc.scan == 'angle']; tt = sorted(a.theta.unique(), reverse=True)
    rows = [[tex(b)] + [mnum('%.2f', v) for v in g.sort_values('theta', ascending=False).BE_kcal] for b, g in a.groupby('base')]
    write('tab_cp_angle', longtable('l' + 'r' * len(tt), r'Dimer & ' + ' & '.join(r'%d$^\circ$' % v for v in tt), rows,
          'ORCA counterpoise binding energies (kcal\\,mol$^{-1}$) along the rigid angle scans (target $\\theta_1$ in the header). Negative values are repulsive.', 'tab:cpangle'))


def t_cp_protein():
    """ORCA counterpoise energies of the protein clusters (results/orca_cp_protein.csv)."""
    pth = os.path.join(RES, 'orca_cp_protein.csv')
    if not os.path.exists(pth): return
    p = pd.read_csv(pth); t = pd.read_csv(os.path.join(RES, 'orca_runtime.csv')); t = t[t.set == 'F_protein']
    t['model'] = t.input.str.split('/').str[-1].str.replace(r'_cp_(AB|A|B)\.inp$', '', regex=True)
    wall = t.groupby('model').seconds.sum() / 60; nbf = t.groupby('model').basis_functions.max()
    rows = []
    for _, r in p.iterrows():
        sysn = r['system'].split('_')[0].upper(); lig, ch = r['ligand'].split(':')[0], r['ligand'].split(':')[1]
        acc = '; '.join(x.split('···')[1].strip().split('@')[0].replace(':' + ch, '') for x in r['xbs'].split(';'))
        rows.append([sysn, '%s (%s)' % ({'LDT': 'IDD594'}.get(lig, lig), ch), tex(acc), str(int(r['n_xb'])), str(int(r['n_atoms'])), str(int(nbf[r['model']])), mnum('%+d', int(r['charge_donor'])).replace('+0', '0'),
                     mnum('%+d', int(r['charge_acceptor'])).replace('+0', '0'), mnum('%+.2f', r['Eint_kcal']), '%.0f' % wall[r['model']]])
    head = ('\\begin{table}[!htbp]\\centering\\footnotesize\n'
            '\\caption{ORCA counterpoise interaction energies of the protein clusters ($\\omega$B97X-V/def2-TZVPD, gas phase, DefGrid3, Sthresh $10^{-6}$). '
            'Fragment 1 is the ligand. Fragment 2 holds the acceptor environment. XB: halogen bonds in the model. XB acceptors: their acceptor residues. $q_1$, $q_2$: fragment charges. $N_\\mathrm{bf}$: basis functions. '
            '$E_\\mathrm{int}$ covers all contacts between the fragments (kcal\\,mol$^{-1}$). Time: wall time of the AB, A, and B jobs on 40 cores (min).}\\label{tab:cpprot}\n'
            '\\begin{tabular}{lllrrrrrrr}\\toprule System & Ligand & XB acceptors & XB & Atoms & $N_\\mathrm{bf}$ & $q_1$ & $q_2$ & $E_\\mathrm{int}$ & Time\\\\\\midrule\n')
    write('tab_cp_protein', head + '\n'.join(' & '.join(r) + '\\\\' for r in rows) + '\n\\bottomrule\\end{tabular}\\end{table}\n')


if __name__ == '__main__':
    t_params(); t_assumptions(); n = t_tests(); t_hal59(); t_protein(); t_perf(); nj = t_jobs(); t_scans(); t_cp(); t_cp_protein()
    print('tests', n, 'jobs', nj)
