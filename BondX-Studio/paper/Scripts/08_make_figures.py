#!/usr/bin/env python3
"""08_make_figures.py
All data figures of Manuscript_V2 and ESI_V2 (PDF + PNG). The flowchart (Fig. 2) is TikZ: Figures_V2_MS/Fig2_BondX_flowchart.tex.
Figures that need QM outputs (results/qm_nbo.csv, results/qm_eda.csv) are drawn when those files exist;
otherwise a clearly marked placeholder is written so that the LaTeX documents still compile.
Usage: python3 08_make_figures.py ..   (project root)
"""
import os, sys, csv, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
RES = os.path.join(ROOT, 'results'); MS = os.path.join(ROOT, 'Figures_V2_MS'); ESI = os.path.join(ROOT, 'Figures_V2_ESI')
SHOTS = os.path.join(ROOT, 'Figures_V2_MS', 'shots')
for d in (MS, ESI): os.makedirs(d, exist_ok=True)

# reference categorical palette (validated; see dataviz palette): blue, orange, aqua, yellow, magenta, green, violet, red
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
HAL = {'I': (C[0], 'o', 'I'), 'BR': (C[1], 's', 'Br'), 'CL': (C[2], '^', 'Cl'), 'F': ('#8a8a86', 'D', 'F')}
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.edgecolor': '#9a9993', 'axes.linewidth': 0.7, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.titlesize': 9.5, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left',
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True, 'legend.frameon': False,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 300, 'pdf.fonttype': 42})


def save(fig, folder, name):
    fig.savefig(os.path.join(folder, name + '.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(folder, name + '.png'), bbox_inches='tight', dpi=220); plt.close(fig)
    print('wrote', name)


import bxstyle
from bxstyle import save
bxstyle.apply()


def placeholder(folder, name, title, lines, size=(7.0, 2.6)):
    fig, ax = plt.subplots(figsize=size); ax.axis('off')
    ax.add_patch(plt.Rectangle((0.01, 0.03), 0.98, 0.94, transform=ax.transAxes, fill=False, ls='--', lw=1.2, ec='#e34948'))
    ax.text(0.5, 0.78, title, ha='center', va='center', fontsize=11.5, weight='bold', color='#e34948', transform=ax.transAxes)
    ax.text(0.5, 0.40, '\n'.join(lines), ha='center', va='center', fontsize=11.1, color=INK2, transform=ax.transAxes)
    save(fig, folder, name)


# ----------------------------------------------------------------------------- Fig 1: literature
def fig1():
    yrs = list(range(2003, 2027))
    allxb = [11, 16, 21, 19, 30, 48, 58, 61, 95, 116, 153, 161, 200, 188, 255, 222, 255, 245, 281, 257, 270, 232, 290, 280]
    prot = [5, 3, 9, 5, 7, 11, 11, 12, 27, 33, 30, 27, 42, 46, 55, 65, 70, 62, 59, 56, 54, 54, 86, 85]
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ax.bar(yrs[:-1], allxb[:-1], color=C[0], width=0.78, label='All halogen-bond papers')
    ax.bar(yrs[-1], allxb[-1], color=C[0], alpha=0.4, width=0.78, label='2026 (to 18 Nov)')
    ax.plot(yrs, prot, color=C[1], lw=1.8, marker='o', ms=3.2, mec='white', mew=0.6, label='With protein, drug or ligand')
    ax.set_xlabel('Publication year'); ax.set_ylabel('Papers per year'); ax.set_ylim(0, 420); ax.grid(axis='x', visible=False)
    ax.legend(loc='upper left', handlelength=1.4, borderaxespad=0.2)
    save(fig, MS, 'Fig1_literature')


# ----------------------------------------------------------------------------- Fig 3: interface montage
def fig3():
    from PIL import Image
    sd = os.path.join(MS, 'shots')
    def main_area(f, h_css):
        im = Image.open(os.path.join(sd, f)).convert('RGB'); sc = im.width / 1440
        return im.crop((int(262 * sc), int(62 * sc), im.width - int(12 * sc), int(h_css * sc)))
    A = main_area('b_detect.png', 640); Bm = main_area('c_qm_setup.png', 660)
    Cm = Image.open(os.path.join(sd, 'd_qm_result.png')).convert('RGB'); Cm = Cm.crop((0, 0, Cm.width, min(Cm.height, int(Cm.width * 1.25))))
    ra, rb = A.height / A.width, Bm.height / Bm.width
    fig = plt.figure(figsize=(7.2, 4.65 * (ra + rb) + 0.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], height_ratios=[ra, rb], wspace=0.04, hspace=0.12)
    for (r, c, im, t) in ((0, 0, A, '(a) Step 2: geometric detection'), (1, 0, Bm, '(b) Step 4: QM/MM set-up'),
                          (slice(0, 2), 1, Cm, '(c) Step 4: model and charges')):
        ax = fig.add_subplot(gs[r, c]); ax.imshow(im); ax.axis('off'); ax.set_title(t, fontsize=10, loc='left', fontweight='bold')
    save(fig, MS, 'Fig3_interface')


def esi_shots():
    from PIL import Image
    files = sorted(f for f in os.listdir(SHOTS) if f[:2].isdigit())
    titles = {'01': 'Step 1: load structure', '02': 'Step 2: geometry (T1)', '03': 'Step 3: dynamics (T2)', '04': 'Step 4: QM cluster model', '05': 'Step 5: NBO gate',
              '06': 'Step 6: ALMO-EDA', '07': 'Step 7: score and report', '08': 'All parameters', '09': 'Method page', '10': 'Step 4: QM/MM model (1US0)', '11': 'Tutorial'}
    for part, chunk in enumerate([files[:6], files[6:]]):
        rows = (len(chunk) + 1) // 2
        fig, axs = plt.subplots(rows, 2, figsize=(7.2, 2.72 * rows))
        for ax in np.ravel(axs): ax.axis('off')
        for ax, f in zip(np.ravel(axs), chunk):
            im = Image.open(os.path.join(SHOTS, f)).convert('RGB'); w, h = im.size
            ax.imshow(im.crop((0, 0, w, min(h, int(w * 0.72))))); ax.set_title(titles.get(f[:2], f), fontsize=10, fontweight='bold', loc='left')
        plt.subplots_adjust(wspace=0.04, hspace=0.14)
        save(fig, ESI, 'FigS10_screens' if part == 0 else 'Extra_screens_2')


# ----------------------------------------------------------------------------- Fig 4: T1 benchmark
def fig4():
    d = pd.read_csv(os.path.join(RES, 't1_hal59.csv')); s = pd.read_csv(os.path.join(RES, 't1_sensitivity.csv'))
    xb = d[d.role == 'halogen bond']
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.5), gridspec_kw={'width_ratios': [1.2, 1]})
    acc_m = {'N': 'o', 'O': 's', 'S': 'D', 'pi': 'h'}
    for x, (col, _, lab) in HAL.items():
        if x == 'F': continue
        q = xb[(xb.halogen == x) & xb.covered]
        for acc, m in acc_m.items():
            r = q[q.acc_el == acc]
            a.scatter(r.ref_kcal, r.G, s=26, marker=m, color=col, edgecolor='white', linewidth=0.6, zorder=3)
    nc = xb[~xb.covered]
    a.scatter(nc.ref_kcal, np.full(len(nc), -0.06), s=22, marker='x', color='#8a8a86', linewidth=1.0, zorder=3)
    a.set_xscale('log'); a.set_xticks([0.5, 1, 2, 5, 10, 20]); a.set_xticklabels(['0.5', '1', '2', '5', '10', '20']); a.minorticks_off()
    a.set_xlabel('CCSD(T)/CBS binding energy (kcal/mol)'); a.set_ylabel('BondX geometric score G'); a.set_ylim(-0.12, 1.05)
    a.set_title('(a) HAL59 dimers')
    from scipy.stats import spearmanr
    cv = xb[xb.covered]; rho = spearmanr(cv.G, cv.ref_kcal).correlation
    a.text(0.03, 0.97, 'ρ = %.2f (n = %d)' % (rho, len(cv)), transform=a.transAxes, ha='left', va='top', color=INK)
    h1 = [Line2D([], [], ls='', marker='o', color=HAL[k][0], label=HAL[k][2]) for k in ('I', 'BR', 'CL')]
    h2 = [Line2D([], [], ls='', marker=m, mfc='white', mec=INK2, label=('π' if k == 'pi' else k)) for k, m in acc_m.items()]
    a.legend(handles=h1 + h2, ncol=2, loc='lower right', bbox_to_anchor=(1.0, 0.1), columnspacing=0.8, handletextpad=0.3)
    piv = s.pivot(index='theta1Min', columns='cutoff_shift', values='protein_XB'); lit = s.pivot(index='theta1Min', columns='cutoff_shift', values='literature_found')
    im = b.imshow(piv.values, cmap='Blues', aspect='auto', origin='lower', vmin=0, vmax=piv.values.max() * 1.15)
    b.set_xticks(range(piv.shape[1])); b.set_xticklabels(['%+.1f' % v for v in piv.columns]); b.set_yticks(range(piv.shape[0])); b.set_yticklabels(piv.index)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]; b.text(j, i, '%d\n%d/4' % (v, lit.values[i, j]), ha='center', va='center', fontsize=9, linespacing=1.05, color='white' if v > piv.values.max() * 0.6 else INK)
    i0 = list(piv.index).index(140); j0 = list(piv.columns).index(0.0)
    b.add_patch(plt.Rectangle((j0 - 0.5, i0 - 0.5), 1, 1, fill=False, ec=C[1], lw=1.8))
    b.set_xlabel('Shift of all distance cut-offs (Å)'); b.set_ylabel('θ₁ minimum (°)'); b.grid(False)
    b.set_title('(b) Protein test set')
    plt.tight_layout(w_pad=1.5)
    save(fig, MS, 'Fig5_T1_benchmark')


# ----------------------------------------------------------------------------- Fig 5: embedding ESP
SYS_LABEL = {'1us0_idd594': '1US0 (Br···OG1 Thr113)', '2w0s': '2W0S (Br···OD1 Asn65)', '1p5e': '1P5E (Br11···O Leu83)', '4agl': '4AGL (I···O Leu145)',
             '4gql': '4GQL (Cl···O Pro232)', '4dst': '4DST (Cl···N Gly75)'}
MAIN_LIG = {'1us0_idd594': 'LDT:A:320', '2w0s': 'BVP:A:1207', '1p5e': 'TBS:A:301', '4agl': 'P84:A:400', '4gql': 'R47:A:306', '4dst': '9LI:A:201'}
VAR = [('qmmm_ff14SB', 'ff14SB'), ('qmmm_charmm36', 'CHARMM36'), ('qmmm_Z1', 'Z1'), ('qmmm_Z2', 'Z2'), ('qmmm_Z3', 'Z3'),
       ('qmmm_r10', '10 Å'), ('qmmm_r20', '20 Å'), ('qmmm_r25', '25 Å'), ('qmmm_full', 'all')]


def fig5(field=False):
    d = pd.read_csv(os.path.join(RES, 'protein_models.csv'))
    col = 'field_AX_VperA' if field else 'esp_X_V'
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 5.4), sharex=True)
    groups = [C[0]] + [C[0]] + [C[1]] * 3 + [C[2]] * 4
    for ax, (sysn, lab) in zip(np.ravel(axs), SYS_LABEL.items()):
        r = d[(d.system == sysn) & (d.ligand == MAIN_LIG[sysn])].set_index('variant')
        vals = [r.loc[v, col] if v in r.index else np.nan for v, _ in VAR]
        bars = ax.bar(range(len(VAR)), vals, color=groups, width=0.72)
        bars[0].set_edgecolor(INK); bars[0].set_linewidth(1.0)
        ax.axhline(0, color='#9a9993', lw=0.7); ax.set_title(lab, fontsize=10.5)
        ax.axhline(vals[0], color=INK, lw=0.7, ls=':')
        ax.set_xticks(range(len(VAR))); ax.set_xticklabels([l for _, l in VAR], rotation=60, fontsize=9.5, ha='right'); ax.grid(axis='x', visible=False)
    for ax in axs[:, 0]: ax.set_ylabel('Field along A→X (V/Å)' if field else 'MM potential at X (V)')
    hs = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (C[0], C[1], C[2])]
    fig.legend(hs, ['force field (RC, 15 Å)', 'boundary scheme (ff14SB, 15 Å)', 'embedding radius (ff14SB, RC)'], loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.04), columnspacing=1.0, handlelength=1.2)
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, MS if field else ESI, 'Fig6_embedding_field' if field else 'FigS4_embedding_ESP')


# ----------------------------------------------------------------------------- ESI figures
def figS1():
    import subprocess
    P = json.loads(subprocess.run(['node', '-e', "console.log(JSON.stringify(require('%s').defaults()))" % os.path.join(ROOT, 'Scripts', 'engine.js')], capture_output=True, text=True).stdout)
    g, sc = P['geom'], P['score']
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.0)); axs = axs.ravel()
    R = np.linspace(0.7, 1.15, 300)
    for k, x in enumerate(('I', 'BR', 'CL')):
        s = sc['vdw'][x] + sc['vdw']['O']; Rc = g['cutoffs'][x]['O'] / s
        gd = np.clip((Rc - R) / (Rc - sc['Rbest']), 0, 1); axs[0].plot(R, gd, color=HAL[x][0], lw=1.6, label='%s···O, R$_{cut}$ %.3f' % (HAL[x][2], Rc))
    axs[0].set_xlabel('R = d / Σr$_{vdW}$'); axs[0].set_ylabel('term value'); axs[0].legend(loc='upper right'); axs[0].set_title('(a) g$_d$', fontsize=11.1)
    t = np.linspace(120, 180, 200); axs[1].plot(t, np.clip((t - g['theta1Min']) / (180 - g['theta1Min']), 0, 1), color=C[0], lw=1.6)
    axs[1].set_xlabel('θ₁ (C–X···Y, °)'); axs[1].set_title('(b) g$_θ$', fontsize=11.1)
    t2 = np.linspace(60, 180, 200); axs[2].plot(t2, np.exp(-((t2 - sc['th2Opt']) / sc['th2Width']) ** 2), color=C[1], lw=1.6)
    axs[2].set_xlabel('θ₂ (X···Y–R, °)'); axs[2].set_title('(c) g$_2$, lone-pair Y', fontsize=11.1)
    a = np.linspace(0, 90, 200); axs[3].plot(a, np.clip(1 - a / g['piNormalMax'], 0, 1), color=C[2], lw=1.6)
    axs[3].set_xlabel('α (X···centroid vs. ring normal, °)'); axs[3].set_title('(d) g$_2$, π ring', fontsize=11.1)
    for ax in axs: ax.set_ylim(-0.05, 1.08)
    axs[2].set_ylabel('term value'); plt.tight_layout(h_pad=1.2, w_pad=1.5); save(fig, ESI, 'FigS1_score_functions')


def figS2():
    d = pd.read_csv(os.path.join(RES, 't1_scans.csv'))
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.8)); axs = axs.ravel()
    dist = d[d.scan == 'distance']
    groups = [('I', ['21_benI-acetone', '24_benI-Nm3', '26_benI-mSH', '28_CH3I-ben', 'NH3_PhI', 'MeI_pyr', 'FI_pyr']), ('BR', ['20_benBr-acetone', '23_benBr-Nm3', 'NH3_FBr']), ('CL', ['19_benCl-acetone', '22_benCl-Nm3'])]
    for ax, (x, names) in zip(axs[:3], groups):
        for k, n in enumerate(names):
            r = dist[dist.base == n].sort_values('s'); ax.plot(r.s, r.G, marker='o', ms=2.8, lw=1.3, color=C[k % 8], label=n.split('_', 1)[-1] if n[0].isdigit() else n.replace('_', '···'))
        ax.set_title('(%s) %s donors' % ('abc'[['I', 'BR', 'CL'].index(x)], HAL[x][2])); ax.set_xlabel('d / d$_{eq}$'); ax.set_ylabel('G'); ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), handlelength=1.2, fontsize=9)
    axs[0].set_ylabel('G')
    ang = d[d.scan == 'angle']
    for k, n in enumerate(sorted(ang.base.unique())):
        r = ang[ang.base == n].sort_values('theta1'); axs[3].plot(r.theta1, r.G, marker='o', ms=2.8, lw=1.3, color=C[k], label=n.split('_')[1])
    axs[3].axvline(140, color=INK2, ls=':', lw=0.8); axs[3].set_xlabel('θ₁ (°)'); axs[3].set_title('(d) Angle scans'); axs[3].set_ylabel('G'); axs[3].legend(loc='center left', bbox_to_anchor=(1.0, 0.5), handlelength=1.2, fontsize=9)
    plt.tight_layout(w_pad=0.8, h_pad=1.2); save(fig, ESI, 'FigS2_scans_T1')


def figS3():
    d = pd.read_csv(os.path.join(RES, 'protein_t1.csv')); d = d[(d.includeF == False) & (d.status != 'none')]
    cnt = d.groupby(['pdb', 'status']).size().unstack(fill_value=0)
    cnt = cnt.reindex(columns=['XB', 'near-miss'], fill_value=0).sort_values('XB')
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    y = np.arange(len(cnt)); ax.barh(y, cnt['XB'], color=C[0], height=0.62, label='halogen bond (T1)')
    ax.barh(y, cnt['near-miss'], left=cnt['XB'], color='#bdbcb6', height=0.62, label='near-miss')
    ax.set_yticks(y); ax.set_yticklabels([p.upper().replace('_IDD594', '+IDD594').replace('_PROTONATED', ' (polar H)') for p in cnt.index])
    ax.set_xlabel('Contacts (all chains)'); ax.legend(loc='lower right'); ax.grid(axis='y', visible=False)
    plt.tight_layout(); save(fig, ESI, 'Extra_protein_T1')


def figS4():
    d = pd.read_csv(os.path.join(RES, 'performance.csv')); t = d[d.test == 'T1 detect']
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    tiled = t.system.str.contains(' x'); ax.scatter(t[~tiled].atoms, t[~tiled].ms, s=16, color=C[0], edgecolor='white', lw=0.5, label='test-set structures')
    ax.scatter(t[tiled].atoms, t[tiled].ms, s=26, marker='s', color=C[1], edgecolor='white', lw=0.5, label='1P5E tiled ×2, ×4, ×8')
    xs = np.array([500, 1e5]); k = np.median(t.ms / t.atoms); ax.plot(xs, k * xs, color=INK2, ls=':', lw=0.9, label='linear scaling')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Atoms'); ax.set_ylabel('Parse + T1 detection (ms)'); ax.legend(loc='upper left')
    plt.tight_layout(); save(fig, ESI, 'FigS9_performance')


def figS5():
    d = pd.read_csv(os.path.join(RES, 'protein_models.csv'))
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    for k, (sysn, lab) in enumerate(SYS_LABEL.items()):
        r = d[(d.system == sysn) & (d.ligand == MAIN_LIG[sysn])].set_index('variant'); sh = [0, 3, 4, 5]
        n = [r.loc['cluster_pcm4', 'n_atoms']] + [r.loc['shell%d' % s, 'n_atoms'] for s in sh[1:]]
        ax.plot(sh, n, marker='o', ms=3.5, lw=1.5, color=C[k], label=lab.split(' ')[0])
    ax.set_xlabel('Extra shell around X (Å)'); ax.set_ylabel('QM atoms (with H and link H)'); ax.legend(ncol=2, loc='upper left'); ax.set_xticks([0, 3, 4, 5]); ax.set_ylim(0, 290)
    plt.tight_layout(); save(fig, ESI, 'FigS3_qm_size')


def qm_figures():
    nbo = os.path.join(RES, 'qm_nbo.csv'); eda = os.path.join(RES, 'qm_eda.csv')
    if os.path.exists(os.path.join(RES, 'qm_summary.json')):
        print('Fig8, Fig9, FigS7 and FigS8 are drawn by 17_qm_nbo_eda_analysis.py'); return
    if not os.path.exists(nbo):
        placeholder(MS, 'Fig8_NBO_benchmark', 'Figure 8 is drawn from the QM outputs', [
            '(a) NBO ΣE(2) [n(Y)→σ*(C–X)] vs. CCSD(T)/CBS binding energy for the 59 HAL59 dimers (Q-Chem and Gaussian).',
            '(b) ΣE(2) along the distance scans, with the 0.5 kcal/mol gate.', '(c) ALMO-EDA terms and χ_CT for the HAL59 halogen bonds.',
            'Run run_qchem.sh and run_gaussian.sh in qm_package, then Scripts/07_collect_results.py and Scripts/08_make_figures.py.'])
        placeholder(MS, 'Fig9_protein_QM', 'Figure 9 is drawn from the QM outputs', [
            'Protein set: ΣE(2) and gate decision for each contact in cluster, PCM and QM/MM models,',
            'its change with force field, boundary scheme, embedding radius and QM-region size, and the BXQI with tier provenance.'])
        for nm, t in (('FigS5_angle_scan_QM', 'NBO and EDA along the angle scans'), ('FigS6_method_basis', 'Method, basis and program dependence of ΣE(2)')):
            placeholder(ESI, nm, t + ' (from QM outputs)', ['Generated by Scripts/08_make_figures.py after Scripts/07_collect_results.py.'])
        return
    q = pd.read_csv(nbo); hal = pd.read_csv(os.path.join(RES, 't1_hal59.csv'))
    a = q[(q.set == 'A_HAL59')].merge(hal[['name', 'ref_kcal', 'halogen', 'acc_el', 'role']], on='name')
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.5))
    for prog, m in (('qchem', 'o'), ('gaussian', 's')):
        r = a[a.program == prog]
        axs[0].scatter(r.ref_kcal, r.E2_core.clip(lower=0.01), s=16, marker=m, c=[HAL.get(x, ('#888',))[0] for x in r.halogen], edgecolor='white', lw=0.4, label=prog)
    axs[0].set_xscale('log'); axs[0].set_yscale('log'); axs[0].axhline(0.5, color=C[7], ls='--', lw=0.9); axs[0].set_xlabel('CCSD(T)/CBS (kcal/mol)'); axs[0].set_ylabel('ΣE(2) (kcal/mol)'); axs[0].legend(fontsize=8.5)
    axs[0].set_title('(a) HAL59')
    b = q[(q.set == 'B_scan_distance') & (q.program == 'qchem')].copy(); sc = pd.read_csv(os.path.join(RES, 't1_scans.csv')); sc['name'] = sc.file.str[:-4]
    b = b.merge(sc[['name', 'base', 's']], on='name')
    for k, n in enumerate(sorted(b.base.unique())):
        r = b[b.base == n].sort_values('s'); axs[1].plot(r.s, r.E2_core.clip(lower=0.01), lw=1.1, marker='o', ms=2, color=C[k % 8])
    axs[1].set_yscale('log'); axs[1].axhline(0.5, color=C[7], ls='--', lw=0.9); axs[1].set_xlabel('d / d$_{eq}$'); axs[1].set_title('(b) distance scans (Q-Chem)')
    if os.path.exists(eda):
        e = pd.read_csv(eda); e = e[e.set == 'A_HAL59'].merge(hal[['name', 'ref_kcal', 'halogen']], on='name')
        axs[2].scatter(e.ref_kcal * 4.184, -e.total, s=16, c=[HAL.get(x, ('#888',))[0] for x in e.halogen], edgecolor='white', lw=0.4)
        lim = [0, max((e.ref_kcal * 4.184).max(), (-e.total).max()) * 1.05]; axs[2].plot(lim, lim, color=INK2, ls=':', lw=0.8)
        axs[2].set_xlabel('CCSD(T)/CBS (kJ/mol)'); axs[2].set_ylabel('−ΔE$_{int}$ ALMO (kJ/mol)'); axs[2].set_title('(c) ALMO-EDA total')
    plt.tight_layout(); save(fig, MS, 'Fig8_NBO_benchmark')


# ----------------------------------------------------------------------------- Fig 7: ORCA counterpoise benchmark
ACC_M = {'N': 'o', 'O': 's', 'S': 'D', 'pi': 'h', 'P': 'v', 'F': '<'}
def fig_cp():
    pth = os.path.join(RES, 'orca_cp_hal59.csv')
    if not os.path.exists(pth):
        placeholder(MS, 'Fig7_CP_benchmark', 'Figure 7 is drawn from the ORCA counterpoise outputs', ['Run Scripts/07_collect_results.py and Scripts/13_orca_cp_analysis.py.']); return
    h = pd.read_csv(pth); sm = json.load(open(os.path.join(RES, 'orca_cp_summary.json'))); sc = pd.read_csv(os.path.join(RES, 'orca_cp_scans.csv'))
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.2, 3.1), gridspec_kw={'width_ratios': [1.05, 1, 1]})
    lim = [0.18, 30]; a.plot(lim, lim, color=INK2, lw=0.8, ls=':', zorder=1)
    for (X, acc), r in h.groupby(['halogen', 'acc']):
        ctl = r.role != 'halogen bond'; col = HAL[X.upper()][0]
        a.scatter(r.ref_kcal[~ctl], r.BE_kcal[~ctl], s=22, marker=ACC_M[acc], color=col, edgecolor='white', linewidth=0.5, zorder=3)
        a.scatter(r.ref_kcal[ctl], r.BE_kcal[ctl], s=26, marker=ACC_M[acc], facecolor='white', edgecolor=INK2, linewidth=0.9, zorder=3)
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(lim); a.set_ylim(lim); a.minorticks_off()
    for ax_ in (a.xaxis, a.yaxis): ax_.set_ticks([0.3, 1, 3, 10, 30]); ax_.set_ticklabels(['0.3', '1', '3', '10', '30'])
    a.set_xlabel('CCSD(T)/CBS (kcal/mol)'); a.set_ylabel('ORCA CP (kcal/mol)'); a.set_title('(a) HAL59')
    st = sm['hal59_all']
    a.text(0.03, 0.97, 'MAE %.2f kcal/mol\nr = %.3f' % (st['MAE'], st['pearson_r']), transform=a.transAxes, ha='left', va='top', color=INK, linespacing=1.2)
    hh = [Line2D([], [], ls='', marker='o', color=HAL[k][0], label=HAL[k][2]) for k in ('I', 'BR', 'CL')] + [Line2D([], [], ls='', marker='o', mfc='white', mec=INK2, label='control')]
    a.legend(handles=hh, loc='lower right', handletextpad=0.2, borderaxespad=0.1)
    # (b) distance scans, normalised to s = 1.00; open symbols: G = 0 (outside the T1 window)
    d = sc[sc.scan == 'distance']
    for base, r in d.groupby('base'):
        r = r.sort_values('s'); col = HAL[r.halogen.iloc[0].upper()][0]
        b.plot(r.s, r.BE_rel, color=col, lw=1.0, alpha=0.9, zorder=2)
        on = r.G > 0; b.scatter(r.s[on], r.BE_rel[on], s=9, color=col, zorder=3, edgecolor='none')
        b.scatter(r.s[~on], r.BE_rel[~on], s=11, facecolor='white', edgecolor=col, linewidth=0.8, zorder=3)
    b.axhline(0, color=INK2, lw=0.6); b.set_xlabel('d / d$_{eq}$'); b.set_ylabel('BE / BE(d$_{eq}$)'); b.set_ylim(-0.05, 1.08); b.set_title('(b) Distance scans')
    b.legend(handles=[Line2D([], [], ls='', marker='o', ms=4, color=INK2, label='G > 0'), Line2D([], [], ls='', marker='o', ms=4, mfc='white', mec=INK2, label='G = 0')], loc='upper right', handletextpad=0.2)
    # (c) angle scans
    g = sc[sc.scan == 'angle']
    c.axvspan(118, 140, color='#f1efe9', lw=0, zorder=0); c.axhline(0, color=INK2, lw=0.6)
    for base, r in g.groupby('base'):
        r = r.sort_values('theta'); col = HAL[r.halogen.iloc[0].upper()][0]; ls = '-' if r.acc.iloc[0] == 'O' else '--'
        c.plot(r.theta, r.BE_rel, color=col, lw=1.1, ls=ls, marker='o', ms=2.4, zorder=2)
    c.set_xlim(118, 182); c.set_xticks([120, 140, 160, 180]); c.set_xlabel('θ$_1$ (°)'); c.set_ylabel('BE / BE(180°)'); c.set_title('(c) Angle scans')
    c.text(129, 0.75, 'G = 0', ha='center', va='center', color=INK2); c.set_ylim(-1.15, 1.62)
    c.legend(handles=[Line2D([], [], color=INK2, lw=1.1, ls='-', label='O (acetone)'), Line2D([], [], color=INK2, lw=1.1, ls='--', label='N (NMe$_3$)')], loc='upper right', handlelength=1.6, borderaxespad=0.0)
    plt.tight_layout(w_pad=1.0); save(fig, MS, 'Fig7_CP_benchmark')
    # ESI: signed error per dimer
    hs = h.sort_values(['halogen', 'ref_kcal']).reset_index(drop=True)
    nm = lambda n: n.split('_', 1)[1] if n[0].isdigit() else n.replace('_', '···')
    fig, ax = plt.subplots(figsize=(6.4, 10.5))
    y = np.arange(len(hs))[::-1]
    ax.barh(y, hs.err_kcal, color=[HAL[x.upper()][0] if r == 'halogen bond' else '#b9b7b0' for x, r in zip(hs.halogen, hs.role)], height=0.75)
    ax.axvline(0, color=INK2, lw=0.8); ax.set_yticks(y); ax.set_yticklabels([nm(n) for n in hs.name], fontsize=9); ax.set_ylim(-0.8, len(hs) - 0.2)
    ax.set_xlabel('BE(ORCA) − BE(CCSD(T)/CBS) (kcal/mol)'); ax.grid(axis='y', visible=False)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=HAL[k][0], label=HAL[k][2]) for k in ('CL', 'BR', 'I')] + [plt.Rectangle((0, 0), 1, 1, color='#b9b7b0', label='control')], loc='lower left')
    save(fig, ESI, 'Extra_cp_errors')



if __name__ == '__main__':
    ALL = {'fig1': fig1, 'fig3': fig3, 'fig4': fig4, 'fig5': fig5, 'fig6': lambda: fig5(field=True), 'figS1': figS1, 'figS2': figS2, 'figS3': figS3,
           'figS4': figS4, 'figS5': figS5, 'shots': esi_shots, 'cp': fig_cp, 'qm': qm_figures}
    sel = [a for a in sys.argv[2:] if a in ALL] or list(ALL)
    for k in sel: ALL[k]()
