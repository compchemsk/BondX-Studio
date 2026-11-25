#!/usr/bin/env python3
"""17_qm_nbo_eda_analysis.py
Analysis of the Gaussian 16 NBO, Gaussian counterpoise and Q-Chem ALMO-EDA outputs (read by 07_collect_results.py).
Inputs : results/qm_nbo.csv, qm_eda.csv, qm_cp.csv, qm_status.csv, t1_hal59.csv, t1_scans.csv, orca_cp_hal59.csv, orca_cp_scans.csv,
         orca_cp_protein.csv, protein_t1.csv
Outputs: results/qm_hal59.csv, qm_scans.csv, qm_method_basis.csv, qm_eda_settings.csv, qm_protein.csv, qm_summary.json
         Figures_V2_MS/Fig8_NBO_benchmark, Fig9_protein_QM; Figures_V2_ESI/FigS7_angle_scan_QM, FigS8_method_basis
         tables/tab_qm_hal59.tex, tab_qm_eda.tex, tab_qm_prot.tex
Usage: python3 17_qm_nbo_eda_analysis.py ..
Energies: NBO E(2) in kcal/mol; EDA terms in kJ/mol in the tables (1 kcal = 4.184 kJ).
"""
import os, sys, json, math
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bxutil import cli

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
RES = os.path.join(ROOT, 'results'); MS = os.path.join(ROOT, 'Figures_V2_MS'); ESI = os.path.join(ROOT, 'Figures_V2_ESI'); TAB = os.path.join(ROOT, 'tables')
KJ = 4.184; GATE = 0.5
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
HALC = {'CL': C[2], 'BR': C[1], 'I': C[0], 'F': '#8a8a86'}
ACCM = {'N': 'o', 'O': 's', 'S': 'D', 'pi': '^', 'P': 'v', 'F': '<'}
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


def r3(x): return None if x is None or (isinstance(x, float) and math.isnan(x)) else round(float(x), 3)
def mnum(fmt, v):
    s = fmt % v
    return s.replace('-', '$-$') if s.startswith('-') else s

nbo = pd.read_csv(os.path.join(RES, 'qm_nbo.csv')); eda = pd.read_csv(os.path.join(RES, 'qm_eda.csv')); cp = pd.read_csv(os.path.join(RES, 'qm_cp.csv'))
st = pd.read_csv(os.path.join(RES, 'qm_status.csv')); hal = pd.read_csv(os.path.join(RES, 't1_hal59.csv'))
sc = pd.read_csv(os.path.join(RES, 't1_scans.csv')); sc['name'] = sc.file.str[:-4]
ocp = pd.read_csv(os.path.join(RES, 'orca_cp_hal59.csv')); ocs = pd.read_csv(os.path.join(RES, 'orca_cp_scans.csv'))
for d in (nbo, eda, cp): assert d.finished.all()
S = {}

# ------------------------------------------------------------------ job status
gq = st[~st.program.str.startswith('orca')].copy(); gq['ok'] = gq.status == 'ok'
S['status'] = {f'{s}|{p}': {'n': int(len(g)), 'ok': int(g.ok.sum())} for (s, p), g in gq.groupby(['set', 'program'])}
qn = gq[gq.program == 'qchem']; S['qchem_nbo_no_E2_table'] = int(qn.status.str.contains('SECOND ORDER').sum())

# ------------------------------------------------------------------ HAL59
a = nbo[nbo.set == 'A_HAL59'].merge(hal[['name', 'subset', 'ref_kcal', 'role', 'halogen', 'acc_el', 'covered', 'pass_T1', 'G']], on='name', how='left')
assert len(a) == 59 and a.name.is_unique
e = eda[eda.set == 'A_HAL59'].set_index('name'); g = cp[(cp.set == 'A_HAL59') & (cp.program == 'gaussian')].set_index('name')
a['eda_total_kJ'] = a.name.map(e.total); a['eda_frz_kJ'] = a.name.map(e.frz); a['eda_pol_kJ'] = a.name.map(e.pol); a['eda_ct_kJ'] = a.name.map(e.ct)
a['eda_ep_kJ'] = a.name.map(e.elec_pauli); a['eda_disp_kJ'] = a.name.map(e.disp); a['chi_CT'] = a.name.map(e.chi_CT); a['mechanism'] = a.name.map(e.mechanism)
a['g16_cp_kcal'] = a.name.map(g.Eint_CP_kJ) / KJ; a['orca_cp_kcal'] = a.name.map(ocp.set_index('name').Eint_kcal)
a['eda_total_kcal'] = a.eda_total_kJ / KJ
a['control'] = a.role.str.startswith('control')
assert a[['eda_total_kJ', 'g16_cp_kcal', 'orca_cp_kcal']].notna().all().all()
# BXQI from G, E and C (dimers: the EDA total is the contact energy)
bx = []
for _, r in a.iterrows():
    j = {'G': float(r.G), 'E2': float(r.E2_core), 'eda': {'total': float(r.eda_total_kJ)}}
    q = cli('bxqi', '--json', json.dumps(j)); bx.append((q.get('score'), q.get('grade', '')))
a['BXQI'] = [b[0] for b in bx]; a['BXQI_grade'] = [b[1] for b in bx]
a.to_csv(os.path.join(RES, 'qm_hal59.csv'), index=False)
xb = a[~a.control]; cov = xb[xb.covered]
rho_all = spearmanr(xb.ref_kcal, xb.E2_core); rho_cov = spearmanr(cov.ref_kcal, cov.E2_core)
rho_cov_noion = spearmanr(cov[cov.ionicX != True].ref_kcal, cov[cov.ionicX != True].E2_core)
ct_e2 = spearmanr(-xb.eda_ct_kJ, xb.E2_core)
def err(x, ref):
    d = x - ref; return {'MAE': r3(d.abs().mean()), 'RMSD': r3(np.sqrt((d ** 2).mean())), 'MSE': r3(d.mean()), 'MaxAE': r3(d.abs().max()), 'r': r3(pearsonr(x, ref)[0])}
be_eda = -a.eda_total_kcal; be_g16 = -a.g16_cp_kcal; be_orca = -a.orca_cp_kcal
S['hal59'] = {
    'n': 59, 'n_xb': int(len(xb)), 'n_covered': int(len(cov)), 'controls': a[a.control][['name', 'E2_core', 'verified']].values.tolist(),
    'controls_all_rejected': bool((~a[a.control].verified).all()),
    'covered_all_verified': bool(cov.verified.all()), 'xb_verified': int(xb.verified.sum()), 'xb_rejected': xb[~xb.verified][['name', 'acc_el', 'E2_core', 'ref_kcal']].values.tolist(),
    'rho_E2_ref_covered': r3(rho_cov.statistic), 'rho_E2_ref_all_xb': r3(rho_all.statistic), 'rho_E2_ref_covered_without_ionic': r3(rho_cov_noion.statistic),
    'rho_E2_negCT': r3(ct_e2.statistic), 'rho_G_E2_covered': r3(spearmanr(cov.G, cov.E2_core).statistic),
    'E2_min_covered': r3(cov.E2_core.min()), 'E2_max_covered': r3(cov.E2_core.max()),
    'E2_over_BE_median_covered': r3((cov.E2_core / cov.ref_kcal).median()),
    'ionic': a[a.ionicX == True][['name', 'E2_core', 'lp_sigma', 'lp_ry']].values.tolist(),
    'ry_share_median': r3((xb.lp_ry / xb.E2_core.where(xb.E2_core > 0)).median()), 'ry_share_max': r3((xb.lp_ry / xb.E2_core.where(xb.E2_core > 0)).max()),
    'npa_ct_max': r3(xb.ct_NPA.max()), 'rho_E2_npaCT': r3(spearmanr(xb.E2_core, xb.ct_NPA).statistic),
    'err_eda_vs_ref': err(be_eda, a.ref_kcal), 'err_g16cp_vs_ref': err(be_g16, a.ref_kcal), 'err_orca_vs_ref': err(be_orca, a.ref_kcal),
    'eda_minus_orca_maxabs_kcal': r3((a.eda_total_kcal - a.orca_cp_kcal).abs().max()), 'eda_minus_orca_mean_kcal': r3((a.eda_total_kcal - a.orca_cp_kcal).mean()),
    'chi_ct_range_xb': [r3(xb.chi_CT.min()), r3(xb.chi_CT.max())], 'chi_ct_median_xb': r3(xb.chi_CT.median()),
    'mechanism_counts_xb': xb.mechanism.value_counts().to_dict(),
    'ct_share_of_attraction_median': r3((xb.eda_ct_kJ / (xb.eda_pol_kJ + xb.eda_ct_kJ + xb.eda_disp_kJ)).median()),
    'frz_positive_xb': int((xb.eda_frz_kJ > 0).sum()),
    'halogen_series': {k: xb[xb.name.isin(v)].set_index('name').loc[v, 'E2_core'].round(2).tolist() for k, v in
                       {'acetone Cl/Br/I': ['19_benCl-acetone', '20_benBr-acetone', '21_benI-acetone'], 'NMe3 Cl/Br/I': ['22_benCl-Nm3', '23_benBr-Nm3', '24_benI-Nm3']}.items()},
    'bxqi_covered_min_max': [r3(cov.BXQI.min()), r3(cov.BXQI.max())],
    'grades_covered': cov.BXQI_grade.str[0].value_counts().to_dict(),
}

# ------------------------------------------------------------------ distance and angle scans
b = nbo[nbo.set.isin(['B_scan_distance', 'C_scan_angle'])].merge(sc[['name', 'base', 'scan', 's', 'theta_target', 'G', 'halogen', 'acc_el']], on='name', how='left')
eb = eda[eda.set.isin(['B_scan_distance', 'C_scan_angle'])].set_index('name')
for k in ('total', 'ct', 'pol', 'frz', 'chi_CT'): b['eda_' + k] = b.name.map(eb[k])
assert b[['base', 'eda_total']].notna().all().all()
ocs['key'] = ocs.apply(lambda r: (r.base, 'distance', round(r.s, 2)) if r.scan == 'distance' else (r.base, 'angle', int(r.theta)), axis=1)
ok_ = dict(zip(ocs.key, ocs.BE_kcal))
b['key'] = b.apply(lambda r: (r.base, 'distance', round(r.s, 2)) if r.scan == 'distance' else (r.base, 'angle', int(r.theta_target)), axis=1)
b['orca_BE'] = b.key.map(ok_); assert b.orca_BE.notna().all()
def rel(df, ref_sel, col):
    out = {}
    for base, g_ in df.groupby('base'):
        r0 = g_[ref_sel(g_)][col].iloc[0]; out.update({i: v / r0 for i, v in zip(g_.index, g_[col])})
    return pd.Series(out)
dist = b[b.scan == 'distance'].copy(); ang = b[b.scan == 'angle'].copy()
dist['E2_rel'] = rel(dist, lambda g_: g_.s == 1.0, 'E2_core'); dist['CT_rel'] = rel(dist, lambda g_: g_.s == 1.0, 'eda_ct'); dist['BE_rel'] = rel(dist, lambda g_: g_.s == 1.0, 'orca_BE')
ang['E2_rel'] = rel(ang, lambda g_: g_.theta_target == 180, 'E2_core'); ang['CT_rel'] = rel(ang, lambda g_: g_.theta_target == 180, 'eda_ct'); ang['BE_rel'] = rel(ang, lambda g_: g_.theta_target == 180, 'orca_BE')
pd.concat([dist, ang]).drop(columns=['key']).to_csv(os.path.join(RES, 'qm_scans.csv'), index=False)
def rng(df, sel, col): v = df[sel(df)][col]; return [r3(v.min()), r3(v.max())]
first_fail = {}
for base, g_ in dist.groupby('base'):
    f = g_[(g_.s > 1.0) & (g_.E2_core < GATE)].s; first_fail[base] = r3(f.min()) if len(f) else None
first_fail_ang = {}
for base, g_ in ang.groupby('base'):
    f = g_[g_.E2_core < GATE].theta_target; first_fail_ang[base] = int(f.max()) if len(f) else None
S['scans'] = {
    'dist_E2_rel_1_2': rng(dist, lambda d: d.s == 1.2, 'E2_rel'), 'dist_E2_rel_1_5': rng(dist, lambda d: d.s == 1.5, 'E2_rel'),
    'dist_CT_rel_1_2': rng(dist, lambda d: d.s == 1.2, 'CT_rel'), 'dist_CT_rel_1_5': rng(dist, lambda d: d.s == 1.5, 'CT_rel'),
    'dist_BE_rel_1_2': rng(dist, lambda d: d.s == 1.2, 'BE_rel'), 'dist_E2_rel_0_9': rng(dist, lambda d: d.s == 0.9, 'E2_rel'),
    'dist_E2_monotonic_all': bool(all((g_.sort_values('s').E2_core.diff().dropna() <= 1e-9).all() for _, g_ in dist.groupby('base'))),
    'dist_gate_first_fail': first_fail, 'dist_gate_fail_count_at_1_5': int((dist[dist.s == 1.5].E2_core < GATE).sum()), 'dist_n_dimers': int(dist.base.nunique()),
    'dist_G0_but_verified': int(((dist.G == 0) & dist.verified).sum()), 'dist_G0_points': int((dist.G == 0).sum()),
    'ang_E2_rel_150': rng(ang, lambda d: d.theta_target == 150, 'E2_rel'), 'ang_E2_rel_140': rng(ang, lambda d: d.theta_target == 140, 'E2_rel'),
    'ang_E2_rel_130': rng(ang, lambda d: d.theta_target == 130, 'E2_rel'), 'ang_E2_rel_120': rng(ang, lambda d: d.theta_target == 120, 'E2_rel'),
    'ang_CT_rel_150': rng(ang, lambda d: d.theta_target == 150, 'CT_rel'), 'ang_CT_rel_120': rng(ang, lambda d: d.theta_target == 120, 'CT_rel'),
    'ang_gate_largest_failing_theta': first_fail_ang, 'ang_verified_at_130': int(ang[ang.theta_target == 130].verified.sum()), 'ang_verified_at_120': int(ang[ang.theta_target == 120].verified.sum()),
    'ang_n_dimers': int(ang.base.nunique()), 'ang_E2_monotonic_all': bool(all((g_.sort_values('theta_target').E2_core.diff().dropna() >= -1e-9).all() for _, g_ in ang.groupby('base'))),
    'ang_ct_attractive_at_120': int((ang[ang.theta_target == 120].eda_ct < 0).sum()), 'ang_orca_repulsive_at_130': int((ang[ang.theta_target == 130].orca_BE < 0).sum()),
}

# ------------------------------------------------------------------ method and basis (Gaussian)
d = nbo[nbo.set == 'D_method_basis'].copy(); d['base'] = d.name.str.split('__').str[0]; d['level'] = d.name.str.split('__', n=1).str[1]
d.to_csv(os.path.join(RES, 'qm_method_basis.csv'), index=False)
pv = d.pivot_table(index='base', columns='level', values='E2_core')
fun = ['wB97X-D__def2-TZVP', 'M06-2X__def2-TZVP', 'B3LYP-D3BJ__def2-TZVP', 'PBE0-D3BJ__def2-TZVP']
rel_f = pv[fun].div(pv['wB97X-D__def2-TZVP'], axis=0)
S['method_basis'] = {'n_dimers': int(pv.shape[0]), 'levels': list(pv.columns),
    'dft_rel_range': [r3(rel_f.values.min()), r3(rel_f.values.max())], 'hf_rel_range': [r3((pv['HF__def2-TZVP'] / pv['wB97X-D__def2-TZVP']).min()), r3((pv['HF__def2-TZVP'] / pv['wB97X-D__def2-TZVP']).max())],
    'svp_rel_range': [r3((pv['wB97X-D__def2-SVP'] / pv['wB97X-D__def2-TZVP']).min()), r3((pv['wB97X-D__def2-SVP'] / pv['wB97X-D__def2-TZVP']).max())],
    'qzvp_rel_range': [r3((pv['wB97X-D__def2-QZVP'] / pv['wB97X-D__def2-TZVP']).min()), r3((pv['wB97X-D__def2-QZVP'] / pv['wB97X-D__def2-TZVP']).max())],
    'gate_changes': int((d.groupby('base').verified.nunique() > 1).sum()),
    'qzvp_rel_sigma_only': [r3(x) for x in (lambda q: [q.min(), q.max()])((d.assign(sig=d.lp_sigma + d.pi_sigma).pivot_table(index='base', columns='level', values='sig')).pipe(lambda t: t['wB97X-D__def2-QZVP'] / t['wB97X-D__def2-TZVP']))],
    'svp_rel_sigma_only': [r3(x) for x in (lambda q: [q.min(), q.max()])((d.assign(sig=d.lp_sigma + d.pi_sigma).pivot_table(index='base', columns='level', values='sig')).pipe(lambda t: t['wB97X-D__def2-SVP'] / t['wB97X-D__def2-TZVP']))],
    'ry_qzvp_max': r3(d[d.level == 'wB97X-D__def2-QZVP'].lp_ry.max()), 'ry_tzvp_max_nonionic': r3(d[(d.level == 'wB97X-D__def2-TZVP') & (d.ionicX != True)].lp_ry.max()),
    'rank_rho_min_vs_wb97xd_tzvp': r3(min(spearmanr(pv[c], pv['wB97X-D__def2-TZVP']).statistic for c in pv.columns))}

# ------------------------------------------------------------------ EDA settings
es = eda[eda.set == 'E_eda_settings'].copy(); es['base'] = es.name.str.split('__').str[0]; es['setting'] = es.name.str.split('__').str[-1]
ref = eda[eda.set == 'A_HAL59'].set_index('name')
es['ref_total'] = es.base.map(ref.total); es['ref_ct'] = es.base.map(ref.ct)
assert es.ref_total.notna().all()
es.to_csv(os.path.join(RES, 'qm_eda_settings.csv'), index=False)
S['eda_settings'] = {s_: {'dTotal_kJ_range': [r3((g_.total - g_.ref_total).min()), r3((g_.total - g_.ref_total).max())],
                          'dCT_kJ_range': [r3((g_.ct - g_.ref_ct).min()), r3((g_.ct - g_.ref_ct).max())],
                          'CT_ratio_range': [r3((g_.ct / g_.ref_ct).min()), r3((g_.ct / g_.ref_ct).max())], 'n': int(len(g_))}
                     for s_, g_ in es.groupby('setting')}

# ------------------------------------------------------------------ protein
f = nbo[nbo.set == 'F_protein'].copy()
f['sys'] = f.name.str.split('__').str[0]; f['lig'] = f.name.str.split('__').str[1]; f['variant'] = f.name.str.split('__').str[2]
pt = pd.read_csv(os.path.join(RES, 'protein_t1.csv')); pt = pt[pt.status == 'XB'].drop_duplicates(['pdb', 'donor', 'X', 'acceptor', 'accAtom'])
def gkey(lbl):
    dn, ac = [x.strip() for x in lbl.split('···')]; return dn.split('@')[0], dn.split('@')[1], ac.split('@')[0], ac.split('@')[1]
gmap = {(r.pdb, r.donor, r.X, r.acceptor, r.accAtom): r.G for r in pt.itertuples()}
f['G'] = [gmap.get((s_,) + gkey(l)) for s_, l in zip(f.sys, f.label)]
assert f.G.notna().all(), f[f.G.isna()].label.unique()
f['contact'] = f.label.str.replace(' ··· ', '···')
f.to_csv(os.path.join(RES, 'qm_protein.csv'), index=False)
main = f[f.lig.str.endswith('A')]
pv2 = main.pivot_table(index='contact', columns='variant', values='E2_core')
S['protein'] = {'models': int(f.name.nunique()), 'contacts': int(f.contact.nunique()), 'rows': int(len(f)),
    'default_cluster': {r.contact: [r3(r.E2_core), bool(r.verified)] for r in f[f.variant == 'cluster_pcm4'].itertuples()},
    'range_over_variants': {c: [r3(pv2.loc[c].min()), r3(pv2.loc[c].max()), int(pv2.loc[c].notna().sum())] for c in pv2.index},
    'gate_flips': [c for c in pv2.index if (pv2.loc[c].dropna() >= GATE).nunique() > 1],
    'pcm_minus_gas': {c: r3(pv2.loc[c, 'cluster_pcm4'] - pv2.loc[c, 'cluster_gas']) for c in pv2.index},
    'qmmm_vs_pcm': {c: r3(pv2.loc[c, 'qmmm_ff14SB'] - pv2.loc[c, 'cluster_pcm4']) for c in pv2.index},
    'ff_charmm_minus_ff14': {c: r3(pv2.loc[c, 'qmmm_charmm36'] - pv2.loc[c, 'qmmm_ff14SB']) for c in pv2.index},
    'boundary_range': {c: [r3(pv2.loc[c, ['qmmm_ff14SB', 'qmmm_Z1', 'qmmm_Z2', 'qmmm_Z3']].min()), r3(pv2.loc[c, ['qmmm_ff14SB', 'qmmm_Z1', 'qmmm_Z2', 'qmmm_Z3']].max())] for c in pv2.index},
    'radius_range': {c: [r3(pv2.loc[c, ['qmmm_r10', 'qmmm_ff14SB', 'qmmm_r20', 'qmmm_r25', 'qmmm_full']].min()), r3(pv2.loc[c, ['qmmm_r10', 'qmmm_ff14SB', 'qmmm_r20', 'qmmm_r25', 'qmmm_full']].max())] for c in pv2.index},
    'shell_range': {c: [r3(pv2.loc[c, ['cluster_pcm4', 'shell3', 'shell4', 'shell5']].min()), r3(pv2.loc[c, ['cluster_pcm4', 'shell3', 'shell4', 'shell5']].max())] for c in pv2.index},
    'his110': {k: r3(v) for k, v in pv2.loc[[c for c in pv2.index if 'THR113' in c][0], ['qmmm_HIS110_HID', 'qmmm_ff14SB', 'qmmm_HIS110_HIP']].items()},
}
# replicas (chain B / C) at the default cluster model
rep = f[f.variant == 'cluster_pcm4'].copy(); rep['res_pair'] = rep.contact.str.replace(r':[A-Z]@', '@', regex=True).str.replace(r'\d{3,4}', '', regex=True)
S['protein']['replica_pairs'] = {k: [r3(x) for x in g_.E2_core] for k, g_ in rep.groupby(rep.sys + '|' + rep.contact.str.split('@').str[1].str.split('···').str[0] + '|' + rep.contact.str.split('@').str[-1]) if len(g_) > 1}
# EDA and counterpoise of the protein clusters
pe = eda[eda.set == 'F_protein'].copy(); pe['sys'] = pe.name.str.split('__').str[0]; pe['lig'] = pe.name.str.split('__').str[1]
gc = cp[(cp.set == 'F_protein') & (cp.program == 'gaussian')].set_index('name'); oc = cp[(cp.set == 'F_protein') & (cp.program == 'orca')].set_index('name')
pe['g16_cp_kJ'] = pe.name.str.replace('_eda', '').map(gc.Eint_CP_kJ); pe['orca_cp_kJ'] = pe.name.str.replace('_eda', '').map(oc.Eint_CP_kJ)
pe.to_csv(os.path.join(RES, 'qm_protein_eda.csv'), index=False)
S['protein']['eda'] = {r.name.replace('__cluster_pcm4_eda', ''): {'frz': r3(r.frz), 'pol': r3(r.pol), 'ct': r3(r.ct), 'total': r3(r.total), 'chi': r3(r.chi_CT), 'orca_cp': r3(r.orca_cp_kJ), 'g16_cp': r3(r.g16_cp_kJ)} for r in pe.itertuples()}
S['protein']['eda_minus_orca_kJ'] = r3((pe.total - pe.orca_cp_kJ).abs().max())
S['protein']['g16cp'] = {k.replace('__cluster_pcm4', ''): r3(v) for k, v in gc.Eint_CP_kJ.items()}
# BXQI of the protein contacts (G and E; the cluster EDA covers every contact of the two fragments)
bq = []
for r in f[f.variant == 'cluster_pcm4'].itertuples():
    q = cli('bxqi', '--json', json.dumps({'G': float(r.G), 'E2': float(r.E2_core)})); bq.append({'contact': r.contact, 'sys': r.sys, 'G': r3(r.G), 'E2': r3(r.E2_core), 'BXQI': r3(q['score']), 'grade': q['grade'], 'provenance': q['provenance']})
bq = pd.DataFrame(bq); bq.to_csv(os.path.join(RES, 'qm_protein_bxqi.csv'), index=False)
S['protein']['bxqi'] = bq[['contact', 'BXQI', 'grade']].values.tolist()
# ------------------------------------------------------------------ wall times (from the outputs)
import re, glob
PKG = os.path.join(ROOT, 'qm_package')
def wall(fn):
    t = open(fn, errors='replace').read()
    if fn.endswith('.log'):
        m = re.findall(r'Elapsed time:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+([\d.]+) seconds', t)
        return sum(int(a) * 86400 + int(b) * 3600 + int(c) * 60 + float(dd) for a, b, c, dd in m) if 'Normal termination' in t else None
    m = re.findall(r'Total job time:\s+([\d.]+)s\(wall\)', t); return float(m[-1]) if m and 'Thank you very much' in t else None
tw = {}
for prog, pat in (('gaussian', '*/gaussian/*.log'), ('gaussian_cp', '*/gaussian_cp/*.log'), ('qchem_eda', '*/qchem_eda/*.out')):
    v = [w for w in (wall(x) for x in glob.glob(os.path.join(PKG, pat))) if w]; tw[prog] = {'n': len(v), 'hours': r3(sum(v) / 3600)}
pe_t = {os.path.basename(x).replace('__cluster_pcm4_eda.out', ''): wall(x) for x in glob.glob(os.path.join(PKG, 'F_protein', 'qchem_eda', '*.out'))}
pe_t = {k: v for k, v in pe_t.items() if v}
rej = [k for k, v in S['protein']['eda'].items()]
gate_fail_clusters = sorted({r.name.replace('__cluster_pcm4', '') for r in f[f.variant == 'cluster_pcm4'].itertuples()} - {r.name.replace('__cluster_pcm4', '') for r in f[(f.variant == 'cluster_pcm4') & f.verified].itertuples()})
tw['protein_eda_hours'] = r3(sum(pe_t.values()) / 3600); tw['protein_eda_n'] = len(pe_t)
tw['protein_clusters_all_rejected'] = gate_fail_clusters
tw['protein_eda_hours_rejected_clusters'] = r3(sum(v for k, v in pe_t.items() if k in gate_fail_clusters) / 3600)
tw['protein_gaussian_nbo_hours'] = r3(sum(w for w in (wall(x) for x in glob.glob(os.path.join(PKG, 'F_protein', 'gaussian', '*.log'))) if w) / 3600)
S['timing'] = tw
if tw['gaussian']['n'] == 0 and os.path.exists(os.path.join(RES, 'qm_summary.json')):  # no QM outputs here (e.g. the GitHub copy): keep the stored wall times
    S['timing'] = json.load(open(os.path.join(RES, 'qm_summary.json'))).get('timing', tw)
json.dump(S, open(os.path.join(RES, 'qm_summary.json'), 'w'), indent=1, default=str)

# ================================================================== figures
# Fig 8: HAL59 NBO, scans, EDA
fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.3)); fig.subplots_adjust(wspace=0.55)
ax = axs[0]
for _, r in a.iterrows():
    acc = 'pi' if r.acc_el == 'pi' else r.acc_el; mk = ACCM.get(acc, 'o'); col = HALC.get(r.halogen, MUTE)
    y = max(r.E2_core, 0.012)
    if r.control: ax.plot(r.ref_kcal, y, marker=mk, ms=4.5, mfc='white', mec=INK2, mew=0.9, ls='none')
    else: ax.plot(r.ref_kcal, y, marker=mk, ms=4.3 if acc in ('N', 'O', 'S', 'pi') else 4.0, mfc=col if r.covered else 'white', mec=col, mew=0.9, ls='none')
ax.axhline(GATE, color=C[7], ls='--', lw=0.9); ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(0.01, 150)
ax.set_xlabel('CCSD(T)/CBS (kcal/mol)'); ax.set_ylabel(r'$\Sigma E^{(2)}$ (kcal/mol)'); ax.set_title('(a) HAL59')
ax.text(0.03, 0.97, 'ρ = %.2f' % S['hal59']['rho_E2_ref_covered'], transform=ax.transAxes, va='top')
hh = [Line2D([], [], ls='none', marker='o', mfc=HALC[k], mec=HALC[k], ms=4, label=l) for k, l in (('CL', 'Cl'), ('BR', 'Br'), ('I', 'I'))] + \
     [Line2D([], [], ls='none', marker='o', mfc='white', mec=INK2, ms=4, label='control / no cut-off')]
fig.legend(handles=hh, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.04), handletextpad=0.2)
ax = axs[1]
for base, g_ in dist.groupby('base'):
    g_ = g_.sort_values('s'); col = HALC.get(g_.halogen.iloc[0], MUTE)
    ax.plot(g_.s, g_.E2_core.clip(lower=0.01), color=col, lw=1.1, marker='o', ms=2.4, alpha=0.9)
ax.axhline(GATE, color=C[7], ls='--', lw=0.9); ax.set_yscale('log'); ax.set_xlabel(r'$d/d_\mathrm{eq}$'); ax.set_ylabel(r'$\Sigma E^{(2)}$ (kcal/mol)'); ax.set_title('(b) Distance scans')
ax = axs[2]
for _, r in a.iterrows():
    acc = 'pi' if r.acc_el == 'pi' else r.acc_el; mk = ACCM.get(acc, 'o'); col = HALC.get(r.halogen, MUTE)
    ax.plot(r.ref_kcal, -r.eda_total_kcal, marker=mk, ms=4, mfc='white' if r.control else col, mec=INK2 if r.control else col, mew=0.8, ls='none')
lim = [0, 22]; ax.plot(lim, lim, color=INK2, ls=':', lw=0.8); ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xticks([0, 5, 10, 15, 20]); ax.set_yticks([0, 5, 10, 15, 20])
ax.set_xlabel('CCSD(T)/CBS (kcal/mol)'); ax.set_ylabel(r'$-\Delta E_\mathrm{int}$, EDA (kcal/mol)'); ax.set_title('(c) ALMO-EDA')
ax.text(0.03, 0.97, 'MAE %.2f\nkcal/mol' % S['hal59']['err_eda_vs_ref']['MAE'], transform=ax.transAxes, va='top')
save(fig, MS, 'Fig8_NBO_benchmark')

# Fig 9: protein
VORD = [('cluster_gas', 'gas'), ('cluster_pcm4', 'PCM ε=4'), ('shell3', '+3 Å'), ('shell4', '+4 Å'), ('shell5', '+5 Å'), ('qmmm_ff14SB', 'ff14SB'), ('qmmm_charmm36', 'CHARMM36'),
        ('qmmm_Z1', 'Z1'), ('qmmm_Z2', 'Z2'), ('qmmm_Z3', 'Z3'), ('qmmm_r10', '10 Å'), ('qmmm_r20', '20 Å'), ('qmmm_r25', '25 Å'), ('qmmm_full', 'all MM')]
CORD = ['BVP1207:A@BR···ASN65:A@OD1', 'TBS301:A@BR11···LEU83:A@O', 'P84400:A@I15···LEU145:A@O', 'TBS301:A@BR12···GLU81:A@O', 'R47306:A@CL1···PRO232:A@O',
        'TBS301:A@BR13···PHE80:A@π(ring)', 'LDT320:A@BR8···THR113:A@OG1', 'TBS301:A@BR10···ILE10:A@O', '9LI201:A@CL2···GLY75:A@N']
CLAB = {'BVP1207:A@BR···ASN65:A@OD1': '2W0S Br···OD1 Asn65', 'TBS301:A@BR11···LEU83:A@O': '1P5E Br11···O Leu83', 'P84400:A@I15···LEU145:A@O': '4AGL I···O Leu145',
        'TBS301:A@BR12···GLU81:A@O': '1P5E Br12···O Glu81', 'R47306:A@CL1···PRO232:A@O': '4GQL Cl···O Pro232', 'TBS301:A@BR13···PHE80:A@π(ring)': '1P5E Br13···π Phe80',
        'LDT320:A@BR8···THR113:A@OG1': '1US0 Br···OG1 Thr113', 'TBS301:A@BR10···ILE10:A@O': '1P5E Br10···O Ile10', '9LI201:A@CL2···GLY75:A@N': '4DST Cl···N Gly75'}
assert set(CORD) == set(pv2.index), set(pv2.index) ^ set(CORD)
M = pv2.loc[CORD, [v for v, _ in VORD]].values
fig = plt.figure(figsize=(7.4, 7.6)); gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], width_ratios=[4.6, 1], hspace=0.42, wspace=0.06)
ax = fig.add_subplot(gs[0, 0])
from matplotlib.colors import LinearSegmentedColormap
cm = LinearSegmentedColormap.from_list('b', ['#f4f7fb', '#a9c8ef', C[0], '#163f73'])
ax.imshow(np.log10(np.clip(M, 0.03, None)), cmap=cm, vmin=np.log10(0.03), vmax=np.log10(6), aspect='auto')
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        ax.text(j, i, '%.2f' % v, ha='center', va='center', fontsize=8.6, color='white' if v > 1.5 else INK, weight='bold' if v < GATE else 'normal')
        if v < GATE: ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=C[7], lw=1.2))
ax.set_xticks(range(len(VORD))); ax.set_xticklabels([l for _, l in VORD], fontsize=9.5, rotation=50, ha='right'); ax.set_yticks(range(len(CORD))); ax.set_yticklabels([CLAB[c] for c in CORD], fontsize=9.8)
ax.grid(False); ax.set_title(r'(a) $\Sigma E^{(2)}$ (kcal/mol) per contact and model', pad=22)
for x in (1.5, 4.5, 6.5, 9.5): ax.axvline(x, color='white', lw=2.5)
for xc, t in ((0.5, 'cluster'), (3.0, 'QM size'), (5.5, 'force field'), (8.0, 'boundary'), (11.5, 'MM radius')):
    ax.text(xc, -0.75, t, ha='center', va='bottom', fontsize=9.5, color=INK2)
# (b) BXQI, same rows as (a)
axq = fig.add_subplot(gs[0, 1], sharey=ax)
bb = bq.set_index('contact').loc[CORD]
cols = [C[7] if 'rejected' in g_ else (C[0] if g_[0] in 'AB' else C[3]) for g_ in bb.grade]
axq.barh(range(len(CORD)), bb.BXQI, color=cols, height=0.62)
for yi, (v, g_) in enumerate(zip(bb.BXQI, bb.grade)): axq.text(v + 3, yi, g_.split(' ')[0], va='center', fontsize=9.5)
axq.set_xlim(0, 120); axq.set_xticks([0, 50, 100]); axq.tick_params(axis='y', left=False, labelleft=False); axq.grid(axis='y', visible=False)
axq.set_xlabel('BXQI'); axq.set_title('(b) BXQI', pad=22)
# (c) EDA
ax = fig.add_subplot(gs[1, :])
pe_ = pe.copy(); pe_['lab'] = pe_.sys.str.upper().str.replace('_IDD594', '') + ' ' + pe_.lig.str[:-1] + ' (chain ' + pe_.lig.str[-1] + ')'
pe_ = pe_.sort_values('total'); y = np.arange(len(pe_))
for k, (col, lab) in enumerate((('frz', 'FRZ'), ('pol', 'POL'), ('ct', 'CT'))):
    ax.barh(y + (k - 1) * 0.26, pe_[col], height=0.25, color=[MUTE, C[2], C[0]][k], label=lab)
ax.plot(pe_.total, y, ls='none', marker='D', ms=6, mfc=INK, mec='white', label='total')
ax.axvline(0, color=INK2, lw=0.8); ax.set_yticks(y); ax.set_yticklabels(pe_.lab); ax.grid(axis='y', visible=False)
ax.set_xlabel('ALMO-EDA term (kJ/mol, gas phase)'); ax.set_title('(c) ALMO-EDA of six clusters'); ax.legend(loc='lower right', ncol=4, handlelength=1.2)
ax.set_xlim(-50, 35); ax.set_ylim(-1.3, len(pe_) - 0.5)
save(fig, MS, 'Fig9_protein_QM')

# Fig S7: angle scans (NBO, CT, ORCA)
fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.3), sharey=True); fig.subplots_adjust(wspace=0.12, top=0.80)
for ax, col, lab in zip(axs, ('E2_rel', 'CT_rel', 'BE_rel'), ('Value relative to 180°', None, None)):
    for base, g_ in ang.groupby('base'):
        g_ = g_.sort_values('theta_target', ascending=False); ax.plot(g_.theta_target, g_[col], marker='o', ms=2.6, lw=1.2, color=HALC.get(g_.halogen.iloc[0], MUTE))
    ax.axvspan(118, 140, color='#f3f2ee', zorder=0, lw=0); ax.set_xlim(182, 118); ax.axhline(0, color=INK2, lw=0.6); ax.set_xlabel(r'$\theta_1$ (deg)'); ax.set_xticks([180, 160, 140, 120]); ax.set_ylim(-1.15, 1.1)
    if lab: ax.set_ylabel(lab)
axs[0].set_title(r'(a) NBO $\Sigma E^{(2)}$'); axs[1].set_title('(b) EDA CT'); axs[2].set_title('(c) ORCA binding')
fig.legend(handles=[Line2D([], [], color=HALC[h], marker='o', ms=4, lw=1.5, label=h.title()) for h in HALC if h in set(ang.halogen)], loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.0))
save(fig, ESI, 'FigS7_angle_scan_QM')

# Fig S8: method / basis
LV = [('wB97X-D__def2-SVP', 'ωB97X-D/SVP'), ('wB97X-D__def2-TZVP', 'ωB97X-D/TZVP'), ('wB97X-D__def2-QZVP', 'ωB97X-D/QZVP'), ('M06-2X__def2-TZVP', 'M06-2X'),
      ('B3LYP-D3BJ__def2-TZVP', 'B3LYP-D3(BJ)'), ('PBE0-D3BJ__def2-TZVP', 'PBE0-D3(BJ)'), ('HF__def2-TZVP', 'HF')]
fig, ax = plt.subplots(figsize=(7.2, 4.4))
hx = d.drop_duplicates('base').set_index('base').name
bh = sc.drop_duplicates('base').set_index('base').halogen
order = pv.loc[:, 'wB97X-D__def2-TZVP'].sort_values().index
for k, (lv, lab) in enumerate(LV):
    ax.plot(np.arange(len(order)) + (k - 3) * 0.1, pv.loc[order, lv], ls='none', marker='o^sDvP<'[k], ms=3.8, color=[C[6], C[0], '#163f73', C[1], C[3], C[2], INK2][k], label=lab)
ax.axhline(GATE, color=C[7], ls='--', lw=0.9); ax.set_yscale('log'); ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=50, ha='right', fontsize=9.2)
ax.set_ylabel(r'$\Sigma E^{(2)}$ (kcal/mol)'); ax.legend(ncol=4, loc='lower center', bbox_to_anchor=(0.5, 1.0), columnspacing=1.0, handletextpad=0.2); ax.grid(axis='x', visible=False)
save(fig, ESI, 'FigS8_method_basis')

# ================================================================== tables
def tex(s): return str(s).replace('_', r'\_').replace('···', r'$\cdots$').replace('π', r'$\pi$').replace('#', r'\#')
def lt(spec, head, rows, cap, lab, size=r'\scriptsize'):
    o = [r'{%s' % size, r'\setlength{\tabcolsep}{3pt}', r'\begin{longtable}{%s}' % spec, r'\caption{%s}\label{%s}\\' % (cap, lab), r'\toprule', head + r'\\', r'\midrule', r'\endfirsthead',
         r'\toprule', head + r'\\', r'\midrule', r'\endhead', r'\bottomrule', r'\endlastfoot'] + [' & '.join(r) + r'\\' for r in rows] + [r'\end{longtable}}']
    return '\n'.join(o) + '\n'
rows = []
for _, r in a.iterrows():
    acc = r'$\pi$' if r.acc_el == 'pi' else r.acc_el
    gate = 'pass' if r.verified else 'fail'
    rows.append([tex(r['name']) + (r'$^c$' if r.control else '') + (r'$^i$' if r.ionicX == True else ''), r.halogen.title(), acc, '%.2f' % r.ref_kcal, '%.2f' % r.G, '%.2f' % r.E2_core, '%.2f' % r.lp_ry if r.acc_el != 'pi' else '--',
                 mnum('%.3f', r.ct_NPA), gate, mnum('%.2f', -r.g16_cp_kcal), mnum('%.2f', -r.eda_total_kcal), '%.0f' % r.BXQI if r.BXQI is not None else '--'])
open(os.path.join(TAB, 'tab_qm_hal59.tex'), 'w').write(lt('lllrrrrrlrrr', r'Dimer & X & Y & $E_\mathrm{ref}$ & $G$ & $\Sigma\ET$ & RY & $\Delta q_\mathrm{CT}$ & Gate & BE(G16) & BE(EDA) & BXQI',
    rows, r'NBO gate on the 59 HAL59 dimers. Gaussian~16, $\omega$B97X-D/def2-TZVP, gas phase. $E_\mathrm{ref}$: CCSD(T)/CBS binding energy (kcal\,mol$^{-1}$). $\Sigma\ET$: gate sum (kcal\,mol$^{-1}$) with the Rydberg term RY. $\Delta q_\mathrm{CT}$: NPA charge transfer to fragment 2 (e). Gate: $\Sigma\ET \ge 0.5$. BE(G16): counterpoise-corrected binding energy from Gaussian ($\omega$B97X-D/def2-TZVP). BE(EDA): $-\Delta E_\mathrm{int}$ from ALMO-EDA ($\omega$B97X-V/def2-TZVPD). BXQI from $G$, $\Sigma\ET$ and the EDA energy. $^c$Control. $^i$Ionic NBO Lewis structure: donation into the vacant LP$^*$ of the halogen.', 'tab:qmhal'))
rows = []
for _, r in a.iterrows():
    rows.append([tex(r['name']) + (r'$^c$' if r.control else ''), mnum('%.2f', r.eda_frz_kJ), mnum('%.2f', r.eda_ep_kJ), mnum('%.2f', r.eda_disp_kJ), mnum('%.2f', r.eda_pol_kJ), mnum('%.2f', r.eda_ct_kJ), mnum('%.2f', r.eda_total_kJ),
                 '%.2f' % r.chi_CT if not math.isnan(r.chi_CT) else '--'])
open(os.path.join(TAB, 'tab_qm_eda.tex'), 'w').write(lt('lrrrrrrr', r'Dimer & FRZ & ELEC+PAULI & DISP & POL & CT & Total & $\chi_\mathrm{CT}$',
    rows, r'ALMO-EDA of the 59 HAL59 dimers (Q-Chem~6.2, $\omega$B97X-V/def2-TZVPD, gas phase, kJ\,mol$^{-1}$). FRZ = ELEC+PAULI+DISP. Q-Chem joins ELEC and PAULI for systems with an ECP. CT is counterpoise-corrected. $\chi_\mathrm{CT}=E_\mathrm{CT}/(E_\mathrm{POL}+E_\mathrm{CT})$. $^c$Control.', 'tab:qmeda'))
VV = [('cluster_gas', 'Cluster, gas'), ('cluster_pcm4', 'Cluster, PCM'), ('shell3', 'PCM, +3 Å'), ('shell4', 'PCM, +4 Å'), ('shell5', 'PCM, +5 Å'), ('qmmm_ff14SB', 'QM/MM ff14SB'), ('qmmm_charmm36', 'QM/MM CHARMM36'),
      ('qmmm_Z1', 'QM/MM Z1'), ('qmmm_Z2', 'QM/MM Z2'), ('qmmm_Z3', 'QM/MM Z3'), ('qmmm_r10', 'QM/MM 10 Å'), ('qmmm_r20', 'QM/MM 20 Å'), ('qmmm_r25', 'QM/MM 25 Å'), ('qmmm_full', 'QM/MM all'),
      ('qmmm_HIS110_HID', 'QM/MM His110 HID'), ('qmmm_HIS110_HIP', 'QM/MM His110 HIP')]
rows = []
for c in CORD:
    for v, lab in VV:
        if v in pv2.columns and not math.isnan(pv2.loc[c, v]):
            rr = f[(f.contact == c) & (f.variant == v)].iloc[0]
            rows.append([tex(CLAB[c]).replace('Å', r'\AA{}'), lab.replace('Å', r'\AA{}'), '%.2f' % rr.E2_core, '%.2f' % rr.lp_ry, mnum('%.3f', rr.ct_NPA), 'pass' if rr.verified else 'fail'])
for r in f[~f.lig.str.endswith('A')].itertuples():
    rows.append([tex(r.contact.replace('@', ' ')), 'Cluster, PCM (replica)', '%.2f' % r.E2_core, '%.2f' % r.lp_ry, mnum('%.3f', r.ct_NPA), 'pass' if r.verified else 'fail'])
open(os.path.join(TAB, 'tab_qm_prot.tex'), 'w').write(lt('llrrrl', r'Contact & Model & $\Sigma\ET$ & RY & $\Delta q_\mathrm{CT}$ & Gate', rows,
    r'NBO gate for the protein contacts (Gaussian~16, $\omega$B97X-D/def2-TZVP). PCM: $\varepsilon=4$. QM/MM: RC boundary and \SI{15}{\angstrom} radius by default. The model name gives any other choice. $\Sigma\ET$ and RY in kcal\,mol$^{-1}$. $\Delta q_\mathrm{CT}$: NPA charge transfer between the fragments (e).', 'tab:qmprot'))
print(json.dumps(S, indent=1, default=str)[:6000])
