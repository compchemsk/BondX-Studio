#!/usr/bin/env python3
"""25_v3_toc.py - TOC graphic of Manuscript_V3 (ACS size 3.25 in x 1.75 in), Arial style, drawn from the paper data.
Usage: python3 25_v3_toc.py .."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bxstyle3 import plt, C, VCOL, GREY, INK, INK2, MUTE, save
from matplotlib.patches import Circle, Wedge, FancyArrowPatch, Rectangle
ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
nbo = pd.read_csv(os.path.join(ROOT, 'results', 'v3', 'nbo_v3.csv'))
h = pd.read_csv(os.path.join(ROOT, 'results', 'qm_hal59.csv')).merge(nbo[nbo.set == 'A_HAL59'][['name', 'E2_valence']], on='name')
pc = pd.read_csv(os.path.join(ROOT, 'results', 'v3', 'protein_contacts_table.csv'))
fig = plt.figure(figsize=(3.25, 1.75)); fig.patch.set_facecolor('white')
# molecule sketch -------------------------------------------------------------
a = fig.add_axes([0.0, 0.22, 0.37, 0.56]); a.set_xlim(0, 3.6); a.set_ylim(-1.1, 1.3); a.set_aspect('equal'); a.axis('off')
a.plot([0.25, 0.75], [0.55, 0.0], color='#55595c', lw=1.6, solid_capstyle='round'); a.plot([0.25, 0.75], [-0.55, 0.0], color='#55595c', lw=1.6)
a.add_patch(Circle((0.75, 0), 0.22, fc='#7b7f83', ec=INK, lw=0.6, zorder=3)); a.text(0.75, 0, 'C', color='white', fontsize=6.5, ha='center', va='center', zorder=4, fontweight='bold')
a.plot([0.97, 1.18], [0, 0], color='#55595c', lw=1.8)
a.add_patch(Wedge((1.55, 0), 0.52, -32, 32, width=0.12, fc=C[0], ec='none', zorder=2))
a.add_patch(Wedge((1.55, 0), 0.52, 60, 120, width=0.10, fc='#e8a6a2', ec='none', zorder=2)); a.add_patch(Wedge((1.55, 0), 0.52, 240, 300, width=0.10, fc='#e8a6a2', ec='none', zorder=2))
a.add_patch(Circle((1.55, 0), 0.38, fc='#a34a22', ec=INK, lw=0.6, zorder=3)); a.text(1.55, 0, 'Br', color='white', fontsize=7, ha='center', va='center', zorder=4, fontweight='bold')
a.plot([2.1, 2.78], [0, 0], ls=(0, (1.2, 1.2)), color=INK, lw=1.2)
a.add_patch(Circle((3.0, 0), 0.24, fc='#d94a43', ec=INK, lw=0.6, zorder=3)); a.text(3.0, 0, 'O', color='white', fontsize=6.5, ha='center', va='center', zorder=4, fontweight='bold')
a.plot([3.0, 3.0], [0.25, 0.75], color='#55595c', lw=1.3); a.plot([3.07, 3.07], [0.25, 0.75], color='#55595c', lw=1.3)
a.text(2.2, 0.5, 'σ-hole', color=C[0], fontsize=5.6, ha='center'); a.text(2.5, -0.62, r'$n_\mathrm{O}\!\rightarrow\!\sigma^*_\mathrm{C-Br}$', fontsize=5.6, ha='center', color=INK2)
# data strip ------------------------------------------------------------------
b = fig.add_axes([0.47, 0.36, 0.50, 0.40]); b.set_xscale('log'); b.set_xlim(0.03, 150); b.set_ylim(-0.6, 1.6)
b.axvspan(0.3, 0.8, color=GREY, lw=0, zorder=0)
x = h[(h.role == 'halogen bond')].E2_valence.values; rng = np.random.default_rng(3)
b.scatter(x, 1 + rng.uniform(-0.22, 0.22, len(x)), s=5, color=C[0], alpha=0.8, lw=0, zorder=3)
ctl = h[h.control].E2_valence.clip(lower=0.035).values; b.scatter(ctl, np.full(len(ctl), 1), s=7, marker='x', color=MUTE, lw=0.8, zorder=3)
vc = lambda v: VCOL['supported'] if v >= 0.8 else (VCOL['inconclusive'] if v >= 0.3 else VCOL['not supported'])
pv = pc.val_default.dropna().clip(lower=0.035).values
b.scatter(pv, np.zeros(len(pv)), s=16, marker='D', c=[vc(v) for v in pv], edgecolor='white', lw=0.4, zorder=3)
b.set_yticks([0, 1]); b.set_yticklabels(['protein', 'HAL59'], fontsize=5.8); b.tick_params(axis='x', labelsize=6, pad=1.5, length=2); b.tick_params(axis='y', length=0, pad=2)
b.set_xticks([0.1, 1, 10, 100]); b.set_xticklabels(['0.1', '1', '10', '100'])
b.set_xlabel(r'valence $\Sigma E^{(2)}$ (kcal/mol)', fontsize=5.8, labelpad=1); b.grid(False)
for sp in ['left', 'top', 'right']: b.spines[sp].set_visible(False)
b.text(0.49, 1.52, 'grey zone', fontsize=5.3, ha='center', color='#8a6a00'); b.text(9, 1.52, 'orbital support', fontsize=5.3, ha='center', color=VCOL['supported'])
# tier line -------------------------------------------------------------------
fig.text(0.03, 0.93, 'BondX Studio', fontsize=9, fontweight='bold', color=INK, va='top')
fig.text(0.97, 0.915, 'halogen bonds, screened in the browser', fontsize=5.8, color=INK2, va='top', ha='right')
fig.add_artist(plt.Line2D([0.03, 0.97], [0.82, 0.82], color='#d6d5d0', lw=0.6))
fig.text(0.5, 0.05, 'T1 geometry → T2 MD → T3 QM model, NBO → T4 EDA → T5 rank', fontsize=5.6, ha='center', color=INK)
os.makedirs(os.path.join(ROOT, 'MS_Figures_V3'), exist_ok=True)
fig.savefig(os.path.join(ROOT, 'MS_Figures_V3', 'TOC_graphic.pdf')); fig.savefig(os.path.join(ROOT, 'MS_Figures_V3', 'TOC_graphic.png'), dpi=600); print('wrote TOC_graphic')
