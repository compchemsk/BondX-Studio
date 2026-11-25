#!/usr/bin/env python3
"""16_hiyd_figures_tables.py
Figures and LaTeX tables of the hIYD case study from the CSV/JSON files of 15_hiyd_analysis.py.
  Figures_V2_MS/Fig10_hIYD.pdf          occupancy along the pathways, acceptor classes, NBO re-analysis
  Figures_V2_ESI/FigS10_hiyd_residues.pdf  residues with T1 contacts per pathway (all vs strict)
  Figures_V2_ESI/FigS11_hiyd_nbo.pdf       per-cluster NBO terms (intramolecular vs the T1 pair)
  tables/tab_hiyd_classes.tex, tab_hiyd_nbo.tex, tab_hiyd_residues.tex
Usage: python3 16_hiyd_figures_tables.py ..
"""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
RES = os.path.join(ROOT, 'results', 'hiyd'); MS = os.path.join(ROOT, 'Figures_V2_MS'); ESI = os.path.join(ROOT, 'Figures_V2_ESI'); TAB = os.path.join(ROOT, 'tables')
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
INK, INK2, GRID, MUTE = '#0b0b0b', '#52514e', '#e4e3df', '#9a9993'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.edgecolor': MUTE, 'axes.linewidth': 0.7, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.titlesize': 9.5, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True, 'legend.frameon': False,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 300, 'pdf.fonttype': 42})
def save(fig, folder, name):
    fig.savefig(os.path.join(folder, name + '.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(folder, name + '.png'), bbox_inches='tight', dpi=220); plt.close(fig); print('wrote', name)

import bxstyle
from bxstyle import save
bxstyle.apply()


occ = pd.read_csv(os.path.join(RES, 'hiyd_window_occupancy.csv'))
cls = pd.read_csv(os.path.join(RES, 'hiyd_acceptor_classes.csv'))
nbo = pd.read_csv(os.path.join(RES, 'hiyd_nbo_reanalysis.csv'))
resc = pd.read_csv(os.path.join(RES, 'hiyd_residue_contacts.csv'))
S = json.load(open(os.path.join(RES, 'hiyd_summary.json')))

# ------------------------------------------------------------------ Fig 10
fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.25], hspace=0.48, wspace=0.18, top=0.88)
for j, P in enumerate('ABC'):
    ax = fig.add_subplot(gs[0, j]); d = occ[occ.pathway == P].sort_values('xi')
    ax.axvspan(9.3, 19.3, color='#f3f2ee', zorder=0, lw=0)
    ax.plot(d.xi, 100 * d.pub_frac_10ns, ls='none', marker='o', ms=3.4, mfc='white', mec=MUTE, mew=0.8, label='Published screen (first 10 ns)')
    ax.plot(d.xi, 100 * d.bondx_T1_frac, color=C[0], lw=1.5, marker='o', ms=2.6, mec='white', mew=0.4, label='BondX T1 (25 ns)')
    ax.plot(d.xi, 100 * d.bondx_strict_frac, color=C[1], lw=1.5, marker='s', ms=2.4, mec='white', mew=0.4, label='BondX T1, Lewis bases only')
    ax.set_title(f'({"abc"[j]}) Pathway {P}'); ax.set_xlabel(r'$\xi$ (Å)'); ax.set_xlim(3.3, 19.3); ax.set_ylim(-2, 52)
    ax.set_xticks([4, 8, 12, 16])
    if j == 0: ax.set_ylabel('Frames with contact (%)')
    else: ax.set_yticklabels([])
    if j == 0: fig.legend(*ax.get_legend_handles_labels(), loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.0), handlelength=1.8, columnspacing=1.2)
# (d) acceptor classes
gs2 = gs[1, :].subgridspec(1, 2, width_ratios=[1.45, 1], wspace=0.55)
ax = fig.add_subplot(gs2[0])
order = cls.sort_values('total', ascending=True)
names = {'backbone C=O': 'Backbone C=O', 'backbone N-H': 'Backbone N–H', 'side-chain O': 'Side-chain O', 'aromatic ring': 'Aromatic ring (π)',
         'side-chain N-H (no lone pair)': 'Side-chain N–H', 'side-chain N (lone pair)': 'His N', 'FMN N': 'FMN N5', 'thioether S': 'Met SD'}
y = np.arange(len(order))
ax.barh(y, order.kept_by_strict_filter, color=C[0], height=0.62, label='Lewis base (kept)')
rem = order.total - order.kept_by_strict_filter; m = rem > 0
ax.barh(y[m], rem[m], left=order.kept_by_strict_filter[m], color='white', edgecolor=C[7], hatch='////', lw=0.8, height=0.62, label='No lone pair (flagged weak)')
for yi, t in zip(y, order.total): ax.text(t + 3, yi, str(t), va='center', color=INK)
ax.set_yticks(y); ax.set_yticklabels([names.get(c, c) for c in order['class']]); ax.grid(axis='y', visible=False)
ax.set_xlabel('T1 contacts (11 361 frames)'); ax.set_xlim(0, max(order.total) * 1.18)
ax.set_title('(d) Acceptor classes'); ax.legend(loc='lower right')
# (e) NBO
ax = fig.add_subplot(gs2[1])
col = {'intra': C[1], 'inter': C[0], 'none': MUTE}
for _, r in nbo.iterrows():
    filled = r.lewis == 'strong'
    ax.plot(r.E_pub, r.E2_bondx_pair, ls='none', marker='o' if filled else 'D', ms=4.2 if filled else 3.6,
            mfc=col[r.E_pub_source] if filled else 'white', mec=col[r.E_pub_source], mew=0.9, alpha=0.95)
ax.axhline(0.5, color=INK2, ls=':', lw=0.8); ax.axvline(0.5, color=INK2, ls=':', lw=0.8)
ax.plot([0, 2.8], [0, 2.8], color=GRID, lw=0.8, zorder=0)
ax.set_xlim(-0.1, 2.8); ax.set_ylim(-0.1, 2.8)
ax.set_xlabel('Published NBO (kcal/mol)'); ax.set_ylabel(r'BondX $\Sigma E^{(2)}$ (kcal/mol)')
ax.set_title('(e) NBO, 56 clusters'); ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1, 2])
h = [Line2D([], [], ls='none', marker='o', mfc=C[1], mec=C[1], ms=4, label='Largest: intra'),
     Line2D([], [], ls='none', marker='o', mfc=C[0], mec=C[0], ms=4, label='Largest: inter'),
     Line2D([], [], ls='none', marker='D', mfc='white', mec=INK2, ms=3.6, label='No lone pair')]
ax.legend(handles=h, loc='upper left', bbox_to_anchor=(0.17, 1.0), handletextpad=0.2, borderaxespad=0.2, fontsize=9.2)
save(fig, MS, 'Fig10_hIYD')

# ------------------------------------------------------------------ ESI Fig S10: residues
fig, axs = plt.subplots(1, 3, figsize=(7.2, 4.6))
for ax, P in zip(axs, 'ABC'):
    d = resc[resc.pathway == P].sort_values('frames_T1', ascending=True)
    y = np.arange(len(d))
    ax.barh(y, d.frames_strict, color=C[0], height=0.65)
    rem = (d.frames_T1 - d.frames_strict).values; m = rem > 0
    ax.barh(y[m], rem[m], left=d.frames_strict.values[m], color='white', edgecolor=C[7], hatch='////', lw=0.7, height=0.65)
    ax.set_yticks(y); ax.set_yticklabels(d.residue, fontsize=9.6); ax.grid(axis='y', visible=False); ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_title(f'({"abc"["ABC".index(P)]}) Pathway {P}'); ax.set_xlabel('Frames'); ax.set_xlim(0, 55)
fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.legend(handles=[Patch(color=C[0], label='Lewis-base acceptor'), Patch(facecolor='white', edgecolor=C[7], hatch='////', label='No lone pair')], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.0)); save(fig, ESI, 'FigS10_hiyd_residues')

# ------------------------------------------------------------------ ESI Fig S11: per-cluster terms
d = nbo.copy(); d['lab'] = d.pathway + ': ' + d.residue + ' ' + d.acceptor_atom
fig, axs = plt.subplots(1, 3, figsize=(7.2, 6.0), sharex=True)
for ax, P in zip(axs, 'ABC'):
    e = d[d.pathway == P].reset_index(drop=True); y = np.arange(len(e))[::-1]
    ax.barh(y + 0.19, e.E2_intra_phenolateO, height=0.36, color=C[1], label='LP(phenolate O of I-Tyr) → σ*(C–I)')
    ax.barh(y - 0.19, e.E2_bondx_pair, height=0.36, color=[C[0] if l == 'strong' else 'white' for l in e.lewis], edgecolor=C[0], lw=0.8, label='T1 pair (BondX gate)')
    ax.plot(e.E_pub, y, ls='none', marker='|', ms=7, mew=1.4, color=INK, label='Published value')
    ax.axvline(0.5, color=INK2, ls=':', lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(e.residue + ' ' + e.acceptor_atom, fontsize=9.2); ax.grid(axis='y', visible=False)
    ax.set_title(f'({"abc"["ABC".index(P)]}) Pathway {P}'); ax.set_xlabel(r'$E^{(2)}$ (kcal/mol)')
hh = [Patch(color=C[1], label='LP(phenolate O of I-Tyr) → σ*(C–I)'), Patch(color=C[0], label='T1 pair, Lewis-base acceptor'),
      Patch(facecolor='white', edgecolor=C[0], label='T1 pair, acceptor without lone pair'), Line2D([], [], ls='none', marker='|', ms=7, mew=1.4, color=INK, label='Published value')]
fig.tight_layout(rect=(0, 0, 1, 0.915)); fig.legend(handles=hh, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.0), columnspacing=1.0)
save(fig, ESI, 'FigS11_hiyd_nbo')

# ------------------------------------------------------------------ tables
def tex(s): return str(s).replace('*', '$^{*}$').replace('_', r'\_')
with open(os.path.join(TAB, 'tab_hiyd_classes.tex'), 'w') as f:
    f.write('\\begin{table}[!htbp]\n\\caption{Acceptor classes of the T1 contacts of the iodine of I-Tyr in all 108 umbrella-sampling windows (11\\,361 frames, 25~ns per window). Kept: contacts with a Lewis-base acceptor (strict filter of \\BX{}).}\n\\label{tab:hiydclasses}\n\\centering\\small\n\\begin{tabular}{lrrrrr}\n\\toprule\nAcceptor class & A & B & C & Total & Kept \\\\\n\\midrule\n')
    for _, r in cls.iterrows():
        f.write(f"{names.get(r['class'], r['class']).replace('π', '$\\pi$')} & {r.contacts_A} & {r.contacts_B} & {r.contacts_C} & {r.total} & {r.kept_by_strict_filter} \\\\\n")
    f.write(f"\\midrule\nAll & {cls.contacts_A.sum()} & {cls.contacts_B.sum()} & {cls.contacts_C.sum()} & {cls.total.sum()} & {cls.kept_by_strict_filter.sum()} \\\\\n\\bottomrule\n\\end{{tabular}}\n\\end{{table}}\n")
with open(os.path.join(TAB, 'tab_hiyd_nbo.tex'), 'w') as f:
    f.write('{\\footnotesize\\setlength{\\tabcolsep}{3pt}\n\\begin{longtable}{llrlrrrrrrl}\n')
    f.write('\\caption{NBO re-analysis of the 56 published QM-cluster outputs (B3LYP-D3, CPCM $\\varepsilon=4$). $\\xi$: window (\\AA). $d$ (\\AA) and $\\theta_1$ (deg): T1 geometry of the published pose. Pub.: published value. Intra: LP of the phenolate O of I-Tyr $\\to\\sigma^*$(C--I). Pair: gate value of \\BX{}, LP of the T1 acceptor (or ring bonds) $\\to\\sigma^*$(C--I). Source: donor of the largest single term. L: Lewis class of the acceptor in the protein (s strong, w no lone pair). Energies in kcal\\,mol$^{-1}$.}\\label{tab:hiydnbo}\\\\\n')
    hdr = 'Path & Residue & $\\xi$ & Acc. & $d$ & $\\theta_1$ & Pub. & Intra & Pair & L & Source \\\\\n'
    f.write('\\toprule\n' + hdr + '\\midrule\n\\endfirsthead\n\\toprule\n' + hdr + '\\midrule\n\\endhead\n\\bottomrule\n\\endfoot\n')
    for _, r in nbo.iterrows():
        acc = 'ring' if r.acceptor_type == 'aromatic' else r.acceptor_atom
        pair = f'\\textbf{{{r.E2_bondx_pair:.2f}}}' if r.bondx_verified_lewis else f'{r.E2_bondx_pair:.2f}'
        f.write(f"{r.pathway} & {tex(r.residue)} & {r.xi:.1f} & {acc} & {r.d_A:.2f} & {r.theta1:.1f} & {r.E_pub:.2f} & {r.E2_intra_phenolateO:.2f} & {pair} & {'s' if r.lewis == 'strong' else 'w'} & {r.E_pub_source} \\\\\n")
    f.write('\\end{longtable}}\n')
with open(os.path.join(TAB, 'tab_hiyd_residues.tex'), 'w') as f:
    f.write('\\begin{table}[!htbp]\n\\caption{Residues with T1 contacts to the iodine of I-Tyr. Frames: number of frames with a contact (all T1 contacts / Lewis-base contacts only). $\\xi$: range of windows with a contact (\\AA). $^{*}$Chain B.}\n\\label{tab:hiydres}\n\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}\n')
    f.write('\\begin{tabular}{llrrr|llrrr}\n\\toprule\nPath & Residue & All & Lewis & $\\xi$ & Path & Residue & All & Lewis & $\\xi$ \\\\\n\\midrule\n')
    rows = [r for _, r in resc.iterrows()]; half = (len(rows) + 1) // 2
    for i in range(half):
        a = rows[i]; s = f"{a.pathway} & {tex(a.residue)} & {a.frames_T1} & {a.frames_strict} & {a.xi_min:.1f}--{a.xi_max:.1f}"
        if i + half < len(rows):
            b = rows[i + half]; s += f" & {b.pathway} & {tex(b.residue)} & {b.frames_T1} & {b.frames_strict} & {b.xi_min:.1f}--{b.xi_max:.1f}"
        else: s += ' & & & & &'
        f.write(s + ' \\\\\n')
    f.write('\\bottomrule\n\\end{tabular}\n\\end{table}\n')
print('tables written')
