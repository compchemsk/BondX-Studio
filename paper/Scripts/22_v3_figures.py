#!/usr/bin/env python3
"""22_v3_figures.py - every data figure of Manuscript_V3 and ESI_V3 (Arial style of bxstyle3.py, legends without frames).
Usage: python3 22_v3_figures.py ..      Output: MS_Figures_V3/ and ESI_Figures_V3/ (PDF + PNG).
The workflow chart (Fig. 1) and the TOC graphic are TikZ files; Fig. 2 and Fig. S13 are screenshots (script 23);
Fig. 3 comes from 12_parameter_figure.py."""
import os, sys, re, json, math
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bxstyle3 import plt, C, HALC, VCOL, GREY, INK, INK2, MUTE, save, grey
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
R = os.path.join(ROOT, 'results'); V = os.path.join(R, 'v3')
MS = os.path.join(ROOT, 'MS_Figures_V3'); ESI = os.path.join(ROOT, 'ESI_Figures_V3')
KC = 4.184; LO, HI = 0.3, 0.8
S = json.load(open(os.path.join(V, 'summary_v3.json')))
verdict = lambda e: 'supported' if e >= HI else ('inconclusive' if e >= LO else 'not supported')
SHAPE = {'N': 'o', 'O': 's', 'S': 'D', 'AR': '^', 'π': '^', 'P': 'v', 'F': '<'}
HNAME = {'CL': 'Cl', 'BR': 'Br', 'I': 'I'}
nbo = pd.read_csv(os.path.join(V, 'nbo_v3.csv'))
h = pd.read_csv(os.path.join(R, 'qm_hal59.csv')).merge(nbo[nbo.set == 'A_HAL59'][['name', 'E2_valence', 'E2_RY']], on='name')
t1 = pd.read_csv(os.path.join(R, 't1_hal59.csv'))

def hal_legend(ax, loc='upper left', extra=()):
    hs = [Line2D([], [], marker='o', ls='', color=HALC[x], label=HNAME[x]) for x in ['I', 'BR', 'CL']]
    hs += [Line2D([], [], marker=SHAPE[a], ls='', color=MUTE, mfc='white' if a in 'PF' else MUTE, label=l) for a, l in [('N', 'N'), ('O', 'O'), ('S', 'S'), ('AR', 'ring')]]
    ax.legend(handles=hs + list(extra), loc=loc, ncol=2, handletextpad=0.3, columnspacing=0.8, fontsize=9.5)

# ---------------------------------------------------------------- Fig. 4: T1 on HAL59 and parameter sensitivity
def fig4():
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.4, 4.0), gridspec_kw=dict(width_ratios=[1.15, 1]))
    for r in t1.itertuples():
        if r.role != 'halogen bond': continue
        acc = 'AR' if isinstance(r.ring, str) or (isinstance(r.ring, float) and not math.isnan(r.ring)) else r.acc_el
        if r.covered: a.scatter(r.ref_kcal, r.G, c=HALC[r.halogen], marker=SHAPE.get(acc, 'o'), s=40, edgecolor='white', lw=0.5, zorder=3)
        else: a.scatter(r.ref_kcal, -0.04, c=HALC[r.halogen], marker='x', s=30, lw=1.2, zorder=3)
    a.set_xscale('log'); a.set_xlabel('CCSD(T)/CBS binding energy (kcal/mol)'); a.set_ylabel('geometric score $G$'); a.set_ylim(-0.08, 1.04)
    a.set_title('(a) HAL59 halogen bonds')
    hal_legend(a, 'lower right', [Line2D([], [], marker='x', ls='', color=MUTE, label='P, F acceptor (no cut-off)')])
    a.text(0.03, 0.40, r'Spearman $\rho$ = %.2f [%.2f, %.2f], $n$ = 41' % (S['hal59']['rho_G_cov'], *S['hal59']['rho_G_cov_CI']), transform=a.transAxes, va='top', fontsize=10)
    sens = pd.read_csv(os.path.join(R, 't1_sensitivity.csv'))
    piv = sens.pivot(index='theta1Min', columns='cutoff_shift', values='protein_XB'); rec = sens.pivot(index='theta1Min', columns='cutoff_shift', values='literature_found')
    im = b.imshow(piv.values, cmap='Blues', aspect='auto', origin='lower', vmin=0, vmax=piv.values.max() * 1.25)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            b.text(j, i, f'{piv.values[i, j]}\n{rec.values[i, j]}/4', ha='center', va='center', fontsize=9, color='white' if piv.values[i, j] > 14 else INK)
    b.set_xticks(range(piv.shape[1])); b.set_xticklabels([f'{x:+.1f}' for x in piv.columns]); b.set_yticks(range(piv.shape[0])); b.set_yticklabels([f'{x}°' for x in piv.index])
    j0 = list(piv.columns).index(0.0); i0 = list(piv.index).index(140); b.add_patch(Rectangle((j0 - 0.5, i0 - 0.5), 1, 1, fill=False, ec=C[1], lw=2.2))
    b.set_xlabel('shift of all distance cut-offs (Å)'); b.set_ylabel(r'lower limit of $\theta_1$'); b.grid(False); b.set_title('(b) protein test set: contacts, literature recovered')
    fig.tight_layout(); save(fig, MS, 'Fig4_T1_benchmark')

# ---------------------------------------------------------------- Fig. 5: MM field (explicit H, radius) and its link to SigmaE(2)
SYSN = {'1us0_idd594': '1US0', '2w0s': '2W0S', '1p5e': '1P5E', '4agl': '4AGL', '4gql': '4GQL', '4dst': '4DST'}
def fig5():
    fd = pd.read_csv(os.path.join(V, 'field_v3.csv')); fd['r'] = fd.radius.replace('full', 60).astype(float)
    fig, axs = plt.subplots(1, 3, figsize=(12.2, 3.9), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
    a = axs[0]
    for k, s in enumerate(SYSN):
        d = fd[(fd.system == s) & (fd.ff == 'ff14SB') & (fd.boundary == 'RC')]; d = d[d.contact == d.contact.iloc[0]]
        for eh, ls, mk in [(0, '--', 'o'), (6, '-', 's')]:
            q = d[d.explicitH == eh].sort_values('r'); a.plot(q.r, q.E_AX, ls=ls, marker=mk, ms=4, color=C[k], label=SYSN[s] if eh == 6 else None)
    a.set_xticks([10, 15, 20, 25, 30, 40, 60]); a.set_xticklabels(['10', '15', '20', '25', '30', '40', 'all'])
    a.set_xlabel('embedding radius (Å)'); a.set_ylabel(r'$E_{AX}$ (V Å$^{-1}$)'); a.set_title('(a) field along A→X')
    a.legend(ncol=3, fontsize=9, loc='upper center', bbox_to_anchor=(0.5, -0.2), handlelength=1.6); a.text(0.5, -0.43, 'solid: explicit MM hydrogens; dashed: folded hydrogens', transform=a.transAxes, ha='center', va='top', fontsize=9, color=INK2)
    b = axs[1]; keys = ['explicitH_shift', 'RCD_shift', 'r15_vs_full']; lab = ['explicit H − folded', 'RCD − RC', '15 Å − whole']
    x = np.arange(len(SYSN)); w = 0.26
    for k, (key, l) in enumerate(zip(keys, lab)): b.bar(x + (k - 1) * w, [S['field'][key][s] for s in SYSN], w, color=[C[2], C[6], C[3]][k], label=l)
    b.axhline(0, color=INK2, lw=0.8); b.set_xticks(x); b.set_xticklabels(list(SYSN.values()), rotation=40); b.set_ylabel(r'$\Delta E_{AX}$ (V Å$^{-1}$)'); b.set_title('(b) effect of each setting'); b.legend(fontsize=9, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    c = axs[2]; m = pd.read_csv(os.path.join(V, 'm12_field_vs_dE.csv')); m = m[m['first'] & m.E_AX.notna()]
    lab = {'BVP': '2W0S', 'LDT': '1US0', 'P84': '4AGL', 'R47': '4GQL', 'TBS': '1P5E', '9LI': '4DST'}
    for k, (ct, g) in enumerate(m.groupby('contact')):
        nm = lab[ct[:3]]; col = C[list(SYSN.values()).index(nm)]
        c.scatter(g.E_AX, g.val, color=col, s=26, label=nm, edgecolor='white', lw=0.4)
    c.set_xlabel(r'$E_{AX}$ (V Å$^{-1}$)'); c.set_ylabel(r'valence $\Sigma E^{(2)}$ (kcal/mol)'); c.set_yscale('log'); c.set_title('(c) QM/MM variants'); c.legend(fontsize=9, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.2))
    c.text(0.5, -0.43, r'within a contact: $\rho$ = %.2f–%.2f' % (min(x[1] for x in S['m12']['within'] if x[1] == x[1]), max(x[1] for x in S['m12']['within'] if x[1] == x[1])) + r'; across systems: $\rho$ = %.2f' % S['m12']['rho_across_default'], transform=c.transAxes, va='top', ha='center', fontsize=9, color=INK2)
    fig.tight_layout(); save(fig, MS, 'Fig5_embedding_field')

# ---------------------------------------------------------------- Fig. 6: ORCA counterpoise benchmark
def fig6():
    cp = pd.read_csv(os.path.join(R, 'orca_cp_hal59.csv')); sc = pd.read_csv(os.path.join(R, 'orca_cp_scans.csv'))
    fig, axs = plt.subplots(1, 3, figsize=(12.2, 3.8)); a = axs[0]
    acc = {r.name: ('AR' if isinstance(r.ring, str) else r.acc_el) for r in t1.itertuples()}
    for r in cp.itertuples():
        hal = r.halogen.upper(); ctl = r.role != 'halogen bond'; mk = SHAPE.get(acc.get(r.name, r.acc), 'o')
        a.scatter(r.ref_kcal, r.BE_kcal, marker=mk, s=36, color=HALC.get(hal, MUTE), facecolors='white' if ctl else HALC.get(hal, MUTE), edgecolor=HALC.get(hal, MUTE), lw=0.9, zorder=3)
    lim = [-1, 22]; a.plot(lim, lim, ':', color=INK2, lw=1); a.set_xlim(lim); a.set_ylim(lim)
    a.set_xlabel('CCSD(T)/CBS binding energy (kcal/mol)'); a.set_ylabel('ORCA counterpoise (kcal/mol)'); a.set_title('(a) HAL59')
    a.text(0.04, 0.96, 'MAE 0.40, RMSE 0.54 kcal/mol\nmean signed error −0.07 kcal/mol', transform=a.transAxes, va='top', fontsize=9.5); hal_legend(a, 'lower right')
    for ax, scan, xcol, xl in [(axs[1], 'distance', 's', r'$d/d_\mathrm{eq}$'), (axs[2], 'angle', 'theta', r'$\theta_1$ (°)')]:
        d = sc[sc.scan == scan]
        for k, (bname, g) in enumerate(d.groupby('base')):
            g = g.sort_values(xcol); col = HALC[g.halogen.iloc[0].upper()]
            ax.plot(g[xcol], g.BE_rel, '-', color=col, lw=1.2, alpha=0.9)
            ax.scatter(g[xcol], g.BE_rel, s=18, color=col, facecolors=np.where(g.G > 0, col, 'white'), edgecolor=col, zorder=3)
        ax.axhline(0, color=INK2, lw=0.8); ax.set_xlabel(xl); ax.set_ylabel('binding energy / reference value')
        if scan == 'angle': ax.axvspan(118, 140, color='#eeeeea', zorder=0); ax.invert_xaxis(); ax.set_title(r'(c) angle scans (open: $G$ = 0)')
        else: ax.set_title(r'(b) distance scans (open: $G$ = 0)')
    axs[1].legend(handles=[Line2D([], [], color=HALC[x], label=HNAME[x]) for x in ['I', 'BR', 'CL']], loc='upper right')
    fig.tight_layout(); save(fig, MS, 'Fig6_CP_benchmark')

# ---------------------------------------------------------------- Fig. 7: NBO screen (valence) and ALMO-EDA on the benchmark
def fig7():
    sc = pd.read_csv(os.path.join(R, 'qm_scans.csv')).merge(nbo[nbo.set.isin(['B_scan_distance', 'C_scan_angle'])][['name', 'E2_valence', 'E2_RY']], on='name')
    fig, axs = plt.subplots(1, 3, figsize=(12.4, 4.0)); a = axs[0]
    acc = {r.name: ('AR' if isinstance(r.ring, str) else r.acc_el) for r in t1.itertuples()}
    grey(a)
    for r in h.itertuples():
        hal = r.halogen; open_ = (r.role != 'halogen bond') or (not r.covered); mk = SHAPE.get(acc.get(r.name, r.acc_el), 'o')
        y = max(r.E2_valence, 0.012); col = HALC.get(hal, MUTE)
        a.scatter(max(r.ref_kcal, 0.03), y, marker=mk, s=38, color=col, facecolors='white' if open_ else col, edgecolor=col, lw=0.9, zorder=3)
    a.set_xscale('log'); a.set_yscale('log'); a.set_ylim(0.01, 150); a.set_xlabel('CCSD(T)/CBS binding energy (kcal/mol)'); a.set_ylabel(r'valence $\Sigma E^{(2)}$ (kcal/mol)')
    a.set_title('(a) HAL59'); a.text(0.03, 0.97, r'$\rho$ = %.2f [%.2f, %.2f], $n$ = 41' % (S['hal59']['rho_E_cov'], *S['hal59']['rho_E_cov_CI']), transform=a.transAxes, va='top', fontsize=9.5)
    hal_legend(a, 'lower right', [Line2D([], [], marker='o', ls='', mfc='white', mec=MUTE, color=MUTE, label='open: control, P, or F')])
    b = axs[1]; grey(b, label=False); d = sc[sc.scan == 'distance']
    for bname, g in d.groupby('base'):
        g = g.sort_values('s'); col = HALC[g.halogen.iloc[0]]; b.plot(g.s, g.E2_valence.clip(lower=0.012), '-o', ms=3.5, color=col, lw=1.2)
    b.set_yscale('log'); b.set_ylim(0.01, 150); b.set_xlabel(r'$d/d_\mathrm{eq}$'); b.set_ylabel(r'valence $\Sigma E^{(2)}$ (kcal/mol)'); b.set_title('(b) distance scans')
    b.legend(handles=[Line2D([], [], color=HALC[x], label=HNAME[x]) for x in ['I', 'BR', 'CL']] + [Patch(color=GREY, label='inconclusive 0.3–0.8')], loc='upper right', fontsize=9)
    c = axs[2]; x = h[h.role == 'halogen bond']
    for r in x.itertuples():
        c.scatter(r.ref_kcal, -r.eda_total_kJ / KC, color=HALC[r.halogen], marker=SHAPE.get(acc.get(r.name, r.acc_el), 'o'), s=34, edgecolor='white', lw=0.4, zorder=3)
    lim = [-1, 22]; c.plot(lim, lim, ':', color=INK2, lw=1); c.set_xlim(lim); c.set_ylim(lim); c.set_xlabel('CCSD(T)/CBS binding energy (kcal/mol)'); c.set_ylabel('−ALMO-EDA total (kcal/mol)')
    c.set_title('(c) ALMO-EDA'); c.text(0.04, 0.96, 'MAE 0.39 kcal/mol (59 dimers)', transform=c.transAxes, va='top', fontsize=9.5)
    fig.tight_layout(); save(fig, MS, 'Fig7_NBO_benchmark')

# ---------------------------------------------------------------- Fig. 8: protein contacts
CONTACTS = [('BVP1207:A@BR···ASN65:A@OD1', '2W0S Br···OD1 Asn65'), ('TBS301:A@BR11···LEU83:A@O', '1P5E Br11···O Leu83'), ('P84400:A@I15···LEU145:A@O', '4AGL I···O Leu145'),
            ('TBS301:A@BR12···GLU81:A@O', '1P5E Br12···O Glu81'), ('R47306:A@CL1···PRO232:A@O', '4GQL Cl···O Pro232'), ('TBS301:A@BR13···PHE80:A@π(ring)', '1P5E Br13···π Phe80'),
            ('LDT320:A@BR8···THR113:A@OG1', '1US0 Br···OG1 Thr113*'), ('TBS301:A@BR10···ILE10:A@O', '1P5E Br10···O Ile10'), ('9LI201:A@CL2···GLY75:A@N', '4DST Cl···N Gly75')]
VARS = [('cluster_gas', 'gas'), ('cluster_pcm4', 'PCM ε=4'), ('shell3', '+3 Å'), ('shell4', '+4 Å'), ('shell5', '+5 Å'), ('qmmm_ff14SB', 'ff14SB'), ('qmmm_charmm36', 'CHARMM36'),
        ('qmmm_Z1', 'Z1'), ('qmmm_Z2', 'Z2'), ('qmmm_Z3', 'Z3'), ('qmmm_r10', '10 Å'), ('qmmm_r20', '20 Å'), ('qmmm_r25', '25 Å'), ('qmmm_full', 'all MM')]
def fig8():
    pr = pd.read_csv(os.path.join(R, 'qm_protein.csv')).merge(nbo[nbo.set == 'F_protein'][['name', 'label', 'E2_valence', 'E2_RY', 'hbY', 'hbFrag2']], on=['name', 'label'])
    pr['hbFlag'] = (pr[['hbY', 'hbFrag2']].max(axis=1) >= 0.5) & (pr[['hbY', 'hbFrag2']].max(axis=1) > pr.E2_valence)
    fig = plt.figure(figsize=(12.0, 9.4)); gs = fig.add_gridspec(2, 2, width_ratios=[4.3, 1.25], height_ratios=[1.2, 1], hspace=0.62, wspace=0.07)
    a = fig.add_subplot(gs[0, 0]); M = np.full((len(CONTACTS), len(VARS)), np.nan); FL = set()
    for i, (c, _) in enumerate(CONTACTS):
        for j, (v, _) in enumerate(VARS):
            q = pr[(pr.contact == c) & (pr.variant == v)]
            if len(q): M[i, j] = q.E2_valence.iloc[0]
            if len(q) and q.hbFlag.iloc[0]: FL.add((i, j))
    cmap = plt.get_cmap('Blues'); a.imshow(np.log10(np.clip(M, 0.03, 10)), cmap=cmap, vmin=-1.8, vmax=1.0, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]; vd = verdict(v)
            a.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=9.3, color='white' if v > 1.0 else INK, fontweight='bold' if vd != 'supported' else 'normal')
            if vd != 'supported': a.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=VCOL[vd], lw=1.8))
            if (i, j) in FL: a.text(j + 0.42, i - 0.3, '†', ha='right', va='center', fontsize=10, color='white' if v > 1.0 else INK)
    for xb in [1.5, 4.5, 6.5, 9.5]: a.axvline(xb, color='white', lw=3)
    a.set_xticks(range(len(VARS))); a.set_xticklabels([l for _, l in VARS], rotation=45, ha='right'); a.set_yticks(range(len(CONTACTS))); a.set_yticklabels([l for _, l in CONTACTS]); a.grid(False)
    for x0, l in [(0.5, 'cluster'), (3, 'QM size'), (5.5, 'force field'), (8, 'boundary'), (11.5, 'MM radius')]: a.text(x0, -0.85, l, ha='center', fontsize=10.5, color=INK2)
    a.set_title('(a) valence $\\Sigma E^{(2)}$ (kcal/mol) per contact and model', pad=22)
    a.legend(handles=[Patch(fc='white', ec=VCOL['inconclusive'], lw=1.8, label='inconclusive (0.3–0.8)'), Patch(fc='white', ec=VCOL['not supported'], lw=1.8, label='not supported (< 0.3)'), Line2D([], [], ls='none', marker='$†$', color=INK, ms=9, label='competing H-bond dominates')], loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=3)
    b = fig.add_subplot(gs[0, 1]); pb = pd.read_csv(os.path.join(V, 'bxqi_protein_table.csv'))
    lab = dict(CONTACTS); pb = pb[pb.contact.isin(lab)]; pb['lab'] = pb.contact.map(lab)
    order = [l for _, l in CONTACTS]; pb = pb.set_index('lab').loc[order].reset_index()
    b.barh(range(len(pb)), pb.BXQI, color=[VCOL[v] for v in pb.verdict], height=0.62)
    for i, r in pb.iterrows(): b.text(r.BXQI + 2, i, f'{r.BXQI:.0f}', va='center', fontsize=9.5)
    b.set_yticks(range(len(pb))); b.set_yticklabels([]); b.set_ylim(len(pb) - 0.5, -0.5); b.set_xlim(0, 112); b.set_xlabel('BXQI (G·E)'); b.set_title('(b) ranking index', pad=22)
    b.legend(handles=[Patch(color=VCOL[k], label=k) for k in ['supported', 'inconclusive', 'not supported']], loc='upper center', bbox_to_anchor=(0.45, -0.24), fontsize=9)
    c = fig.add_subplot(gs[1, :]); pe = S['protein_eda']; names = [f"{SYSN[r[0]]} {r[1][:3]} (chain {r[1][3]})" for r in pe]
    y = np.arange(len(pe)); w = 0.26
    for k, (idx, nm, col) in enumerate([(2, 'FRZ', '#9a9993'), (3, 'POL', C[2]), (4, 'CT', C[0])]): c.barh(y + (k - 1) * w, [r[idx] for r in pe], w, color=col, label=nm)
    c.scatter([r[5] for r in pe], y, marker='D', color=INK, s=34, zorder=4, label='total')
    c.axvline(0, color=INK2, lw=0.8); c.set_yticks(y); c.set_yticklabels(names); c.invert_yaxis(); c.set_xlabel('ALMO-EDA term (kcal/mol, gas phase)'); c.legend(ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.2))
    c.set_title('(c) ALMO-EDA of six clusters (every ligand–fragment contact, not the halogen bond alone)')
    save(fig, MS, 'Fig8_protein_QM')

# ---------------------------------------------------------------- Fig. 9: hIYD
def fig9():
    occ = pd.read_csv(os.path.join(R, 'hiyd', 'hiyd_window_occupancy.csv')); cls = pd.read_csv(os.path.join(R, 'hiyd', 'hiyd_acceptor_classes.csv'))
    hn = pd.read_csv(os.path.join(R, 'hiyd', 'hiyd_nbo_reanalysis.csv'))
    fig = plt.figure(figsize=(12.2, 7.6)); gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.32)
    for k, p in enumerate('ABC'):
        a = fig.add_subplot(gs[0, k]); d = occ[occ.pathway == p].sort_values('xi')
        a.axvspan(9.3, 19.3, color='#f2f1ec', zorder=0)
        a.plot(d.xi, d.pub_frac_10ns * 100, 'o', mfc='white', mec=INK2, ms=5, label='published screen (10 ns)')
        a.plot(d.xi, d.bondx_T1_frac * 100, '-', color=C[0], label='BondX T1 (25 ns)'); a.plot(d.xi, d.bondx_strict_frac * 100, '-', color=C[1], label='Lewis-base acceptors')
        a.set_xlabel(r'$\xi$ (Å)'); a.set_ylabel('frames with a contact (%, per window)'); a.set_title(f'({"abc"[k]}) pathway {p}'); a.set_ylim(-2, 45)
        if k == 0: a.legend(loc='upper left', fontsize=9)
    a = fig.add_subplot(gs[1, 0]); cls = cls.sort_values('total'); poor = cls['class'].str.contains('N–H|N-H|no lone|side-chain N-H|backbone N', regex=True)
    cls['class'] = cls['class'].replace({'side-chain N-H (no lone pair)': 'side-chain N–H (poor base)', 'backbone N-H': 'backbone N–H (amide)', 'side-chain N (lone pair)': 'side-chain N (His)'})
    a.barh(cls['class'], cls.total, color=[C[1] if p else C[0] for p in poor], hatch=None)
    for i, (t, p) in enumerate(zip(cls.total, poor)): a.text(t + 3, i, str(t), va='center', fontsize=9.5)
    a.set_xlabel('T1 contacts (all frames)'); a.set_title('(d) acceptor classes'); a.legend(handles=[Patch(color=C[0], label='Lewis base'), Patch(color=C[1], label='poor Lewis base')], loc='lower right', fontsize=9)
    a = fig.add_subplot(gs[1, 1:]); grey(a, 'x')
    src = hn.E_pub_source.fillna('none')
    for s_, col, lab in [('intra', C[1], 'published value = intramolecular O→σ*(C–I)'), ('inter', C[0], 'published value = intermolecular'), ('none', MUTE, 'no term')]:
        q = hn[src == s_]
        for r in q.itertuples(): a.scatter(r.E2_bondx_pair, r.E_pub, marker='D' if r.lewis == 'weak' else 'o', s=36, color=col, facecolors='white' if r.lewis == 'weak' else col, edgecolor=col, lw=1, zorder=3)
        a.scatter([], [], color=col, label=lab)
    a.scatter([], [], marker='D', facecolors='white', edgecolor=INK2, label='poor Lewis base (published model)')
    up = pd.read_csv(os.path.join(V, 'hiyd_decisions_updated.csv')); bb = up[up.source.str.startswith('acetyl')]
    for r in bb.itertuples():
        a.annotate('', xy=(r.E2_used, r.E_pub), xytext=(r.E2_bondx_pair, r.E_pub), arrowprops=dict(arrowstyle='->', color=INK2, lw=0.7, shrinkA=3, shrinkB=2), zorder=2)
        a.scatter(r.E2_used, r.E_pub, marker='D', s=22, color=INK2, zorder=4)
    a.scatter([], [], marker='D', s=22, color=INK2, label='backbone N with the amide restored (acetyl cap)')
    a.set_xlabel(r'BondX valence $\Sigma E^{(2)}$ of the T1 pair (kcal/mol)'); a.set_ylabel('published NBO value (kcal/mol)'); a.set_title('(e) re-analysis of the 56 published NBO outputs')
    a.legend(loc='lower right', fontsize=8.6); a.set_xlim(-0.1, 2.8); a.set_ylim(0.3, 2.8)
    save(fig, MS, 'Fig9_hIYD')

# ---------------------------------------------------------------- ESI figures
def figS1():
    fig, axs = plt.subplots(1, 4, figsize=(12.2, 2.9)); R_ = np.linspace(0.6, 1.2, 300); rc = 3.68 / 3.37
    gd = np.clip((rc - R_) / (rc - 0.85), 0, 1); axs[0].plot(R_, gd, color=C[0]); axs[0].set_xlabel(r'$R = d/(r_X + r_Y)$'); axs[0].set_ylabel(r'$g_d$'); axs[0].set_title('(a) distance (Br···O)')
    t = np.linspace(120, 180, 200); axs[1].plot(t, np.clip((t - 140) / 40, 0, 1), color=C[0]); axs[1].set_xlabel(r'$\theta_1$ (°)'); axs[1].set_ylabel(r'$g_\theta$'); axs[1].set_title('(b) σ-hole angle')
    t2 = np.linspace(60, 180, 200); axs[2].plot(t2, np.exp(-((t2 - 120) / 40) ** 2), color=C[0]); axs[2].set_xlabel(r'$\theta_2$ (°)'); axs[2].set_ylabel(r'$g_2$'); axs[2].set_title('(c) lone-pair acceptor')
    al = np.linspace(0, 90, 200); axs[3].plot(al, np.clip(1 - al / 60, 0, 1), color=C[0]); axs[3].set_xlabel(r'$\alpha$ (°)'); axs[3].set_ylabel(r'$g_2$'); axs[3].set_title('(d) ring acceptor')
    fig.tight_layout(); save(fig, ESI, 'FigS1_score_functions')

def figS2():
    s = pd.read_csv(os.path.join(R, 't1_scans.csv')); fig, axs = plt.subplots(1, 4, figsize=(12.2, 3.0), sharey=True)
    for k, hal in enumerate(['CL', 'BR', 'I']):
        d = s[(s.scan == 'distance') & (s.halogen == hal)]
        for j, (b_, g) in enumerate(d.groupby('base')): g = g.sort_values('s'); axs[k].plot(g.s, g.G, '-o', ms=3.5, color=C[j % 8], label=b_.split('_', 1)[1].replace('Nm3', 'NMe3'))
        axs[k].set_xlabel(r'$d/d_\mathrm{eq}$'); axs[k].set_title(f'({"abc"[k]}) {HNAME[hal]} donors'); axs[k].legend(fontsize=7.5, loc='upper right')
    axs[0].set_ylabel('$G$'); d = s[s.scan == 'angle']
    for j, (b_, g) in enumerate(d.groupby('base')): g = g.sort_values('theta_target'); axs[3].plot(g.theta_target, g.G, '-o', ms=3.5, color=C[j % 8], label=b_.split('_', 1)[1].replace('Nm3', 'NMe3'))
    axs[3].axvline(140, ls=':', color=INK2); axs[3].set_xlabel(r'$\theta_1$ (°)'); axs[3].set_title('(d) angle scans'); axs[3].legend(fontsize=7.5); axs[3].invert_xaxis()
    fig.tight_layout(); save(fig, ESI, 'FigS2_scans_T1')

def figS3_S4():
    pm = pd.read_csv(os.path.join(R, 'protein_models.csv'))
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.0, 3.6), gridspec_kw=dict(width_ratios=[1, 1.5]))
    for k, s in enumerate(SYSN):
        q = pm[(pm.system == s) & pm.variant.isin(['cluster_pcm4', 'shell3', 'shell4', 'shell5'])]
        xs = q.variant.map({'cluster_pcm4': 0, 'shell3': 3, 'shell4': 4, 'shell5': 5}); a.plot(xs, q.n_atoms, '-o', color=C[k], label=SYSN[s])
    a.set_xlabel('extra QM shell (Å)'); a.set_ylabel('QM atoms'); a.set_title('(a) QM-region size'); a.legend(ncol=2, fontsize=9)
    fd = pd.read_csv(os.path.join(V, 'field_v3.csv')); fd['r'] = fd.radius.replace('full', 60).astype(float)
    for k, s in enumerate(SYSN):
        d = fd[(fd.system == s) & (fd.ff == 'ff14SB') & (fd.boundary == 'RC') & (fd.explicitH == 6)]; d = d[d.contact == d.contact.iloc[0]].sort_values('r')
        b.plot(d.r, d.phi_X, '-o', ms=4, color=C[k], label=SYSN[s])
    b.set_xticks([10, 15, 20, 25, 30, 40, 60]); b.set_xticklabels(['10', '15', '20', '25', '30', '40', 'all']); b.set_xlabel('embedding radius (Å)'); b.set_ylabel(r'$\phi_X$ (V)'); b.set_title('(b) MM potential at the halogen'); b.legend(ncol=3, fontsize=9)
    fig.tight_layout(); save(fig, ESI, 'FigS3_size_potential')

def figS4_angles():
    sc = pd.read_csv(os.path.join(R, 'qm_scans.csv')).merge(nbo[nbo.set == 'C_scan_angle'][['name', 'E2_valence', 'E2_RY']], on='name'); d = sc[sc.scan == 'angle']
    fig, axs = plt.subplots(1, 4, figsize=(12.4, 3.1))
    for j, (b_, g) in enumerate(d.groupby('base')):
        g = g.sort_values('theta_target'); r0 = g[g.theta_target == 180].iloc[0]; col = C[j]
        axs[0].plot(g.theta_target, g.E2_valence / r0.E2_valence, '-o', ms=3.5, color=col, label=b_.split('_', 1)[1].replace('Nm3', 'NMe3'))
        axs[1].plot(g.theta_target, 100 * g.E2_RY / (g.E2_RY + g.E2_valence), '-o', ms=3.5, color=col)
        axs[2].plot(g.theta_target, g.eda_ct / r0.eda_ct, '-o', ms=3.5, color=col); axs[3].plot(g.theta_target, g.orca_BE / r0.orca_BE, '-o', ms=3.5, color=col)
    for ax, t, yl in zip(axs, ['(a) valence ΣE(2)', '(b) Rydberg share', '(c) ALMO-EDA CT', '(d) ORCA binding energy'], ['value / value at 180°', 'n→RY* share (%)', 'value / value at 180°', 'value / value at 180°']):
        ax.axvspan(118, 140, color='#eeeeea', zorder=0); ax.invert_xaxis(); ax.set_xlabel(r'$\theta_1$ (°)'); ax.set_title(t); ax.set_ylabel(yl)
    axs[3].axhline(0, color=INK2, lw=0.8); axs[0].legend(fontsize=7.5, loc='lower left'); fig.tight_layout(); save(fig, ESI, 'FigS4_angle_scan_QM')

def figS5_levels():
    mb = pd.read_csv(os.path.join(R, 'qm_method_basis.csv')).merge(nbo[nbo.set == 'D_method_basis'][['name', 'E2_valence']], on='name')
    L = [('wB97X-D__def2-TZVP', 'ωB97X-D/TZVP'), ('M06-2X__def2-TZVP', 'M06-2X/TZVP'), ('B3LYP-D3BJ__def2-TZVP', 'B3LYP-D3(BJ)/TZVP'), ('PBE0-D3BJ__def2-TZVP', 'PBE0-D3(BJ)/TZVP'), ('HF__def2-TZVP', 'HF/TZVP'), ('wB97X-D__def2-SVP', 'ωB97X-D/SVP'), ('wB97X-D__def2-QZVP', 'ωB97X-D/QZVP')]
    bases = sorted(mb.base.unique()); fig, ax = plt.subplots(figsize=(11.5, 3.6)); grey(ax); w = 0.11
    for k, (lv, lab) in enumerate(L):
        q = mb[mb.level == lv].set_index('base').reindex(bases); ax.bar(np.arange(len(bases)) + (k - 3) * w, q.E2_valence, w, color=C[k], label=lab)
    ax.set_yscale('log'); ax.set_xticks(range(len(bases))); ax.set_xticklabels([b.split('_', 1)[1] for b in bases], rotation=35, ha='right'); ax.set_ylabel(r'valence $\Sigma E^{(2)}$ (kcal/mol)')
    ax.legend(ncol=4, fontsize=9, loc='upper left'); ax.set_title('valence ΣE(2) of the 12 scan dimers at seven levels'); fig.tight_layout(); save(fig, ESI, 'FigS5_method_basis')

def figS6_calibration():
    sc = pd.read_csv(os.path.join(R, 'qm_scans.csv')).merge(nbo[nbo.set.isin(['B_scan_distance', 'C_scan_angle'])][['name', 'E2_valence']], on='name')
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.6, 3.9), gridspec_kw=dict(width_ratios=[1.3, 1]))
    grey(a); x = h[h.role == 'halogen bond']; rng_ = np.random.default_rng(1)
    groups = [('HAL59 XB, covered', x[x.covered].E2_valence, C[0]), ('HAL59 XB, P/F', x[~x.covered].E2_valence, C[6]), ('controls', h[h.control].E2_valence, MUTE),
              ('T1-positive scan points,\nCP attractive', sc[(sc.G > 0) & (sc.orca_BE >= 0)].E2_valence, C[2]), ('T1-positive scan points,\nCP repulsive', sc[(sc.G > 0) & (sc.orca_BE < 0)].E2_valence, C[7]),
              ('T1-negative scan points', sc[sc.G == 0].E2_valence, C[3])]
    for k, (lab, v, col) in enumerate(groups): a.scatter(k + rng_.uniform(-0.18, 0.18, len(v)), np.clip(v, 0.012, None), s=16, color=col, alpha=0.85)
    a.set_yscale('log'); a.set_xticks(range(len(groups))); a.set_xticklabels([g[0] for g in groups], rotation=30, ha='right', fontsize=9); a.set_ylabel(r'valence $\Sigma E^{(2)}$ (kcal/mol)')
    a.set_title('(a) decision value in the benchmark sets'); a.legend(loc='upper right', fontsize=9)
    gz = S['greyzone']; rows = [('level of theory (Set D)', gz['level_ratio']), ('protein model variants', gz['model_ratio']), ('combined', gz['combined'])]
    for k, (lab, (lo, hi)) in enumerate(rows): b.plot([lo, hi], [k, k], '-', lw=6, color=C[k], solid_capstyle='butt'); b.text(hi + 0.03, k, f'{lo:.2f}–{hi:.2f}', va='center', fontsize=9.5)
    b.axvline(1, color=INK2, lw=0.8, ls=':'); b.set_yticks(range(3)); b.set_yticklabels([r[0] for r in rows]); b.set_xlim(0.4, 1.9); b.set_ylim(-0.6, 2.6); b.set_xlabel('ΣE(2) relative to the reference model')
    b.set_title('(b) spread used for the grey zone'); b.text(0.42, -0.45, f"0.5/{gz['combined'][1]:.2f} = {gz['zone'][0]:.2f} and 0.5/{gz['combined'][0]:.2f} = {gz['zone'][1]:.2f} → 0.3–0.8 kcal/mol", fontsize=9.5)
    fig.tight_layout(); save(fig, ESI, 'FigS6_gate_calibration')

def figS7_S8_hiyd():
    res = pd.read_csv(os.path.join(R, 'hiyd', 'hiyd_residue_contacts.csv')); fig, axs = plt.subplots(1, 3, figsize=(12.4, 4.2), sharex=False)
    for k, p in enumerate('ABC'):
        d = res[res.pathway == p].sort_values('frames_T1').tail(12); a = axs[k]
        a.barh(d.residue, d.frames_strict, color=C[0], label='Lewis-base contact'); a.barh(d.residue, d.frames_T1 - d.frames_strict, left=d.frames_strict, color='white', ec=C[1], hatch='///', label='poor Lewis base only')
        a.set_title(f'({"abc"[k]}) pathway {p}'); a.set_xlabel('frames with a T1 contact')
    axs[0].legend(fontsize=9, loc='lower right'); fig.tight_layout(); save(fig, ESI, 'FigS7_hiyd_residues')
    hn = pd.read_csv(os.path.join(R, 'hiyd', 'hiyd_nbo_reanalysis.csv')); hn = hn.sort_values(['pathway', 'xi']); fig, ax = plt.subplots(figsize=(12.4, 3.8)); x = np.arange(len(hn))
    ax.axhspan(LO, HI, color=GREY, zorder=0, lw=0)
    ax.bar(x - 0.2, hn.E2_intra_phenolateO, 0.4, color=C[1], label='phenolate O → σ*(C–I) (intramolecular)')
    ax.bar(x + 0.2, hn.E2_bondx_pair, 0.4, color=[C[0] if l == 'strong' else 'white' for l in hn.lewis], ec=C[0], label='BondX valence ΣE(2) of the T1 pair')
    ax.scatter(x, hn.E_pub, marker='_', s=120, color=INK, zorder=4, label='published value')
    up = pd.read_csv(os.path.join(V, 'hiyd_decisions_updated.csv')).sort_values(['pathway', 'xi']).reset_index(drop=True); bbm = up.source.str.startswith('acetyl').values
    ax.scatter(x[bbm] + 0.2, up.E2_used[bbm], marker='D', s=20, color=C[3], zorder=5, label='backbone N, amide restored (acetyl cap)')
    ax.set_xticks(x); ax.set_xticklabels([f'{r.pathway} {r.residue} {r.acceptor_atom}' for r in hn.itertuples()], rotation=90, fontsize=7.5); ax.set_ylabel('E(2) (kcal/mol)'); ax.legend(ncol=2, fontsize=9, loc='upper left')
    ax.set_title('NBO terms of the 56 published cluster outputs (open bars: poor Lewis base)'); ax.set_xlim(-0.8, len(hn) - 0.2); fig.tight_layout(); save(fig, ESI, 'FigS8_hiyd_nbo')

def figS9_t2():
    t2 = pd.read_csv(os.path.join(V, 'hiyd_t2_v3.csv')); fig, axs = plt.subplots(1, 2, figsize=(11.6, 3.7))
    for k, p in enumerate('ABC'):
        d = t2[t2.pathway == p]; axs[0].scatter(d.xi, d.occupancy * 100, s=18 + 60 * d.D / t2.D.max(), color=C[k], alpha=0.8, label=f'pathway {p}', edgecolor='white', lw=0.3)
        axs[1].scatter(d.occupancy * 100, d.tau_ps, s=22, color=C[k], alpha=0.8, label=f'pathway {p}')
    axs[0].axhline(10, ls=':', color=INK2); axs[0].text(18.9, 10.6, 'default report cut-off (10%)', fontsize=9, ha='right'); axs[0].set_xlabel(r'$\xi$ (Å)'); axs[0].set_ylabel('occupancy (%, per window, unweighted)'); axs[0].set_title('(a) T2 occupancy of each contact'); axs[0].legend(fontsize=9)
    for dt in sorted(set(t2.dt_ps)): axs[1].axhline(dt, ls=':', color=MUTE, lw=0.8)
    axs[1].text(axs[1].get_xlim()[1] * 0.98 if False else 25.5, 250, 'frame spacing 160–240 ps', fontsize=8.5, ha='right', color=INK2)
    axs[1].set_xlabel('occupancy (%)'); axs[1].set_ylabel(r'mean residence time $\tau$ (ps)'); axs[1].set_title('(b) residence time'); fig.tight_layout(); save(fig, ESI, 'FigS9_hiyd_T2')

def figS12_classes():
    fig, axs = plt.subplots(1, 2, figsize=(11.6, 5.0), gridspec_kw=dict(width_ratios=[1, 1.25]))
    xb = h[~h.control].copy(); xb['ratio'] = xb.E2_valence / xb.ref_kcal; xb['cls'] = xb.acc_el.replace({'pi': 'ring'})
    order = ['P', 'N', 'S', 'O', 'F', 'ring']; a = axs[0]; rng = np.random.default_rng(3)
    for k, c in enumerate(order):
        d = xb[xb.cls == c]; a.scatter(k + rng.uniform(-0.16, 0.16, len(d)), d.ratio, s=26, color=C[0] if c != 'ring' else C[1], edgecolor='white', lw=0.4, zorder=3)
        a.plot([k - 0.28, k + 0.28], [d.ratio.median()] * 2, color=INK, lw=2, zorder=4)
    a.set_xticks(range(len(order))); a.set_xticklabels(order); a.set_xlabel('acceptor class (HAL59)'); a.set_ylabel(r'valence $\Sigma E^{(2)}$ / binding energy'); a.set_yscale('log'); a.set_yticks([0.4, 0.6, 1, 2, 3, 4]); a.set_yticklabels(['0.4', '0.6', '1', '2', '3', '4']); a.minorticks_off(); a.set_title('(a) orbital donation per unit binding energy')
    ex = pd.read_csv(os.path.join(V, 'extra_nbo.csv')); ex['flags'] = ex['flags'].fillna(''); ex = ex[ex.E2_valence.notna()]
    ex['hb'] = ex[['hbY', 'hbFrag2']].max(axis=1); b = axs[1]
    b.axvspan(LO, HI, color=GREY, zorder=0, lw=0); b.plot([0.01, 6], [0.01, 6], ':', color=INK2, lw=0.8)
    zp = lambda v: np.maximum(v, 0.03)
    G = ex[ex.set == 'G_hard_negatives']; ctrl = G.role.str.startswith('control')
    b.scatter(zp(G[~ctrl].E2_valence), zp(G[~ctrl].hb), marker='D', s=40, facecolors='white', edgecolors=C[1], lw=1.2, label='hIYD frame, poor Lewis base (set G)', zorder=3)
    b.scatter(zp(G[ctrl].E2_valence), zp(G[ctrl].hb), marker='o', s=40, color=C[0], label='hIYD frame, Lewis base (set G)', zorder=3)
    K = ex[ex.set == 'K_hiyd_amide']; b.scatter(zp(K.E2_valence), zp(K.hb), marker='s', s=26, color=MUTE, label='published pose, backbone N with acetyl cap (set K)', zorder=2)
    J = ex[ex.set == 'J_4dst_amide']; b.scatter(zp(J.E2_valence), zp(J.hb), marker='^', s=44, color=C[3], label='4DST Gly75 N, intact amide (set J)', zorder=3)
    b.set_xscale('log'); b.set_yscale('log'); b.set_xlim(0.025, 3); b.set_ylim(0.025, 7)
    b.set_xlabel(r'valence $\Sigma E^{(2)}$ (kcal/mol)'); b.set_ylabel(r'X$\to\sigma^*$(N–H/O–H) (kcal/mol)'); b.set_title('(b) contacts that pass T1 with poor Lewis bases')
    b.text(0.04, 4.6, 'above the diagonal:\ncompeting H-bond dominates', fontsize=8.5, color=INK2, va='top'); b.legend(fontsize=8.5, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    b.set_xticks([0.03, 0.1, 0.3, 1, 3]); b.set_xticklabels(['≤0.03', '0.1', '0.3', '1', '3']); b.set_yticks([0.03, 0.1, 0.3, 1, 3]); b.set_yticklabels(['≤0.03', '0.1', '0.3', '1', '3']); b.minorticks_off()
    fig.tight_layout(); save(fig, ESI, 'FigS12_classes_hardneg')

def figS10_perf():
    p = pd.read_csv(os.path.join(R, 'performance.csv')); d = p[p.test == 'T1 detect']; fig, ax = plt.subplots(figsize=(4.6, 3.4))
    ax.scatter(d.atoms, d.ms, s=20, color=C[0]); xs = np.array([d.atoms.min(), d.atoms.max()]); k = np.median(d.ms / d.atoms); ax.plot(xs, k * xs, ':', color=INK2)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('atoms'); ax.set_ylabel('parse + T1 time (ms)'); fig.tight_layout(); save(fig, ESI, 'FigS10_performance')

if __name__ == '__main__':
    for f in [fig4, fig5, fig6, fig7, fig8, fig9, figS1, figS2, figS3_S4, figS4_angles, figS5_levels, figS6_calibration, figS7_S8_hiyd, figS9_t2, figS10_perf, figS12_classes]: f()
