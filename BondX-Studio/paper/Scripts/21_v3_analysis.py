#!/usr/bin/env python3
"""21_v3_analysis.py - every number of the revised paper (V3) that comes from the existing outputs.
Usage: python3 21_v3_analysis.py ..
Needs results/*.csv of V2 (scripts 03-17) and results/v3/*.csv of scripts 18-20. Writes results/v3/summary_v3.json and
results/v3/*_table.csv. No QM job is run. Spearman confidence intervals: 10000 bootstrap resamples (seed 20261126)."""
import os, sys, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
R = os.path.join(ROOT, 'results'); V = os.path.join(R, 'v3')
KCAL = 4.184
LO, HI = 0.3, 0.8
rng = np.random.default_rng(20261126)
S = {}

def verdict(e): return 'supported' if e >= HI else ('inconclusive' if e >= LO else 'not supported')
def boot_rho(x, y, n=10000):
    x, y = np.asarray(x), np.asarray(y); m = len(x); r = []
    for _ in range(n):
        i = rng.integers(0, m, m); rr = spearmanr(x[i], y[i]).correlation
        if np.isfinite(rr): r.append(rr)
    return [round(float(np.percentile(r, 2.5)), 2), round(float(np.percentile(r, 97.5)), 2)]
def rho(x, y): return round(float(spearmanr(x, y).correlation), 2)
r2 = lambda v: round(float(v), 2)

nbo = pd.read_csv(os.path.join(V, 'nbo_v3.csv'))
# ------------------------------------------------------------------ HAL59
h = pd.read_csv(os.path.join(R, 'qm_hal59.csv'))
hv = nbo[nbo.set == 'A_HAL59'][['name', 'E2_valence', 'E2_RY', 'E2_with_RY', 'hbY']]
h = h.merge(hv, on='name')
h['verdict'] = h.E2_valence.map(verdict)
xb = h[h.role == 'halogen bond']; cov = xb[xb.covered]; pf = xb[~xb.covered]; ctl = h[h.control]
S['hal59'] = dict(n=len(h), n_xb=len(xb), n_cov=len(cov), n_pf=len(pf),
    cov_verdicts=cov.verdict.value_counts().to_dict(), pf_verdicts=pf.verdict.value_counts().to_dict(), xb_verdicts=xb.verdict.value_counts().to_dict(),
    cov_val_range=[r2(cov.E2_valence.min()), r2(cov.E2_valence.max())], cov_inconclusive=cov[cov.verdict != 'supported'][['name', 'E2_valence', 'ref_kcal']].values.tolist(),
    pf_not_supported=pf[pf.verdict != 'supported'][['name', 'E2_valence', 'ref_kcal']].values.tolist(),
    controls=ctl[['name', 'E2_valence', 'verdict']].values.tolist(),
    rho_G_cov=rho(cov.G, cov.ref_kcal), rho_G_cov_CI=boot_rho(cov.G, cov.ref_kcal),
    rho_E_cov=rho(cov.E2_valence, cov.ref_kcal), rho_E_cov_CI=boot_rho(cov.E2_valence, cov.ref_kcal),
    rho_E_all=rho(xb.E2_valence, xb.ref_kcal), rho_E_all_CI=boot_rho(xb.E2_valence, xb.ref_kcal),
    rho_E_NPA=rho(xb.E2_valence, xb.ct_NPA), rho_E_CT=rho(xb.E2_valence, -xb.eda_ct_kJ),
    median_ratio_cov=r2((cov.E2_valence / cov.ref_kcal).median()), median_ratio_all=r2((xb.E2_valence / xb.ref_kcal).median()),
    RY_share_median=r2((xb.E2_RY / xb.E2_with_RY.replace(0, np.nan)).median() * 100), RY_share_max=r2((xb.E2_RY / xb.E2_with_RY.replace(0, np.nan)).max() * 100),
    rho_G_I=rho(cov[cov.halogen == 'I'].G, cov[cov.halogen == 'I'].ref_kcal), rho_G_Br=rho(cov[cov.halogen == 'BR'].G, cov[cov.halogen == 'BR'].ref_kcal))
for acc, nm in [('acetone', 'acetone'), ('Nm3', 'NMe3')]:
    S['hal59']['series_' + nm] = [r2(h[h.name == f'{i}_ben{x}-{acc}'].E2_valence.iloc[0]) for i, x in [(19, 'Cl'), (20, 'Br'), (21, 'I')]] if acc == 'acetone' else [r2(h[h.name == f'{i}_ben{x}-{acc}'].E2_valence.iloc[0]) for i, x in [(22, 'Cl'), (23, 'Br'), (24, 'I')]]
S['hal59']['FI_pyr'] = r2(h[h.name == 'FI_pyr'].E2_valence.iloc[0])
# EDA in kcal/mol
eda = xb.eda_total_kJ / KCAL; orca = xb.orca_cp_kcal; ref = -xb.ref_kcal
S['hal59']['eda_MAE'] = r2((eda - ref).abs().mean()); S['hal59']['eda_vs_orca_max'] = r2((eda - orca).abs().max())
alleda = h.eda_total_kJ / KCAL; S['hal59']['eda_vs_orca_max_all59'] = r2((alleda - h.orca_cp_kcal).abs().max())
S['hal59']['eda_MAE_all59'] = r2((alleda + h.ref_kcal).abs().mean())
S['hal59']['chi_range'] = [r2(xb.chi_CT.min()), r2(xb.chi_CT.max())]; S['hal59']['frz_repulsive'] = int((xb.eda_frz_kJ > 0).sum())
S['hal59']['rho_CT_E'] = rho(-xb.eda_ct_kJ, xb.E2_valence)
# ------------------------------------------------------------------ BXQI (new definition) on HAL59: external check against CCSD(T)/CBS
def bxqi(G, E, C, w=dict(G=0.25, E=0.30, C=0.25), E0=2.0, C0=2.4):
    t = {}
    if G is not None and np.isfinite(G): t['G'] = max(0, min(1, G))
    if E is not None: t['E'] = 1 - np.exp(-max(0, E) / E0)
    if C is not None: t['C'] = 1 - np.exp(-max(0, -C) / C0)
    if any(v <= 0 for v in t.values()): return 0.0
    ws = sum(w[k] for k in t); return float(100 * np.exp(sum(w[k] / ws * np.log(t[k]) for k in t)))
xb = xb.copy()
xb['BXQI'] = [bxqi(r.G if r.covered else None, r.E2_valence, r.eda_total_kJ / KCAL) for r in xb.itertuples()]
xb['BXQI_GE'] = [bxqi(r.G if r.covered else None, r.E2_valence, None) for r in xb.itertuples()]
xb['BXQI_old'] = xb.BXQI_x if 'BXQI_x' in xb else h.set_index('name').loc[xb.name, 'BXQI'].values
S['bxqi_hal59'] = dict(rho_all=rho(xb.BXQI, xb.ref_kcal), rho_all_CI=boot_rho(xb.BXQI, xb.ref_kcal), rho_cov=rho(xb[xb.covered].BXQI, xb[xb.covered].ref_kcal),
    rho_GE_all=rho(xb.BXQI_GE, xb.ref_kcal), rho_old_all=rho(xb.BXQI_old, xb.ref_kcal),
    PF_examples=xb[xb.name.isin(['FI_PCH', 'FI_FMe'])][['name', 'BXQI_old', 'BXQI', 'ref_kcal']].round(1).values.tolist(), range_all=[r2(xb.BXQI.min()), r2(xb.BXQI.max())])
xb[['name', 'ref_kcal', 'covered', 'G', 'E2_valence', 'eda_total_kJ', 'BXQI_old', 'BXQI', 'BXQI_GE']].to_csv(os.path.join(V, 'bxqi_hal59_table.csv'), index=False)
# ------------------------------------------------------------------ scans
sc = pd.read_csv(os.path.join(R, 'qm_scans.csv')).merge(nbo[nbo.set.isin(['B_scan_distance', 'C_scan_angle'])][['name', 'E2_valence', 'E2_RY']], on='name')
dsc = sc[sc.scan == 'distance'].copy(); asc = sc[sc.scan == 'angle'].copy()
ref_d = dsc[dsc.s == 1.0].set_index('base').E2_valence; dsc['rel'] = dsc.E2_valence / dsc.base.map(ref_d)
ref_a = asc[asc.theta_target == 180].set_index('base').E2_valence; asc['rel'] = asc.E2_valence / asc.base.map(ref_a)
S['scans'] = dict(d120_rel=[r2(dsc[dsc.s == 1.2].rel.min() * 100), r2(dsc[dsc.s == 1.2].rel.max() * 100)], d090_rel=[r2(dsc[dsc.s == 0.9].rel.min()), r2(dsc[dsc.s == 0.9].rel.max())],
    d150_verdicts=dsc[dsc.s == 1.5].E2_valence.map(verdict).value_counts().to_dict(), monotone_d=bool(all((g.sort_values('s').E2_valence.diff().dropna() <= 1e-9).all() for _, g in dsc.groupby('base'))),
    a140_rel=[r2(asc[asc.theta_target == 140].rel.min() * 100), r2(asc[asc.theta_target == 140].rel.max() * 100)], a120_rel=[r2(asc[asc.theta_target == 120].rel.min() * 100), r2(asc[asc.theta_target == 120].rel.max() * 100)],
    a120_verdicts=asc[asc.theta_target == 120].E2_valence.map(verdict).value_counts().to_dict(), a130_verdicts=asc[asc.theta_target == 130].E2_valence.map(verdict).value_counts().to_dict(),
    a120_RYshare=[r2((asc[asc.theta_target == 120].E2_RY / (asc[asc.theta_target == 120].E2_RY + asc[asc.theta_target == 120].E2_valence)).min() * 100), r2((asc[asc.theta_target == 120].E2_RY / (asc[asc.theta_target == 120].E2_RY + asc[asc.theta_target == 120].E2_valence)).max() * 100)],
    a180_RYshare_max=r2((asc[asc.theta_target == 180].E2_RY / (asc[asc.theta_target == 180].E2_RY + asc[asc.theta_target == 180].E2_valence)).max() * 100))
# T1-positive benchmark geometries with a repulsive counterpoise energy (specificity against an energy criterion)
t1pos = pd.concat([sc[sc.G > 0][['name', 'scan', 's', 'theta_target', 'G', 'E2_valence', 'orca_BE']]])
rep = t1pos[t1pos.orca_BE < 0]
S['scans']['t1pos_n'] = int(len(t1pos)); S['scans']['t1pos_repulsive'] = rep[['name', 'G', 'E2_valence', 'orca_BE']].round(3).values.tolist()
S['scans']['t1pos_repulsive_verdicts'] = rep.E2_valence.map(verdict).value_counts().to_dict()
S['scans']['t1pos_attractive_verdicts'] = t1pos[t1pos.orca_BE >= 0].E2_valence.map(verdict).value_counts().to_dict()
# ------------------------------------------------------------------ level of theory (Set D) and grey zone
mb = pd.read_csv(os.path.join(R, 'qm_method_basis.csv')).merge(nbo[nbo.set == 'D_method_basis'][['name', 'E2_valence']], on='name')
r0 = mb[mb.level == 'wB97X-D__def2-TZVP'].set_index('base').E2_valence; mb['rel'] = mb.E2_valence / mb.base.map(r0)
lv = mb.groupby('level').rel.agg(['min', 'max']).round(2)
S['levels'] = dict(rel=lv.to_dict(orient='index'), all_range=[r2(mb[mb.level != 'wB97X-D__def2-TZVP'].rel.min()), r2(mb[mb.level != 'wB97X-D__def2-TZVP'].rel.max())],
    verdict_changes=int(sum(len(set(g.E2_valence.map(verdict))) > 1 for _, g in mb.groupby('base'))), min_val=r2(mb.E2_valence.min()),
    rho_min=r2(min(spearmanr(mb[mb.level == l].sort_values('base').E2_valence, mb[mb.level == 'wB97X-D__def2-TZVP'].sort_values('base').E2_valence).correlation for l in mb.level.unique())))
pr = pd.read_csv(os.path.join(R, 'qm_protein.csv')).merge(nbo[nbo.set == 'F_protein'][['name', 'label', 'E2_valence', 'E2_RY', 'E2_with_RY', 'hbY', 'hbFrag2', 'lpOnY']], on=['name', 'label'])
mr = []
for c, d in pr.groupby('contact'):
    if len(d) < 5: continue
    v0 = d[d.variant == 'cluster_pcm4'].E2_valence.iloc[0]
    if v0 < 0.2: continue
    mr += [(d.E2_valence / v0).min(), (d.E2_valence / v0).max()]
lo, hi = min(mr), max(mr); clo, chi = S['levels']['all_range'][0] * lo, S['levels']['all_range'][1] * hi
S['greyzone'] = dict(model_ratio=[r2(lo), r2(hi)], level_ratio=S['levels']['all_range'], combined=[r2(clo), r2(chi)], zone=[r2(0.5 / chi), r2(0.5 / clo)], used=[LO, HI])
# ------------------------------------------------------------------ protein contacts
rows = []
for c, d in pr.groupby('contact'):
    d0 = d[d.variant == 'cluster_pcm4']
    rows.append(dict(contact=c, sys=d.sys.iloc[0], n=len(d), G=r2(d.G.iloc[0]), val_default=r2(d0.E2_valence.iloc[0]) if len(d0) else None, val_min=r2(d.E2_valence.min()), val_max=r2(d.E2_valence.max()),
        RY_default=r2(d0.E2_RY.iloc[0]) if len(d0) else None, withRY_min=r2(d.E2_with_RY.min()), withRY_max=r2(d.E2_with_RY.max()), hbY_max=r2(d.hbY.max()), hbFrag2_max=r2(d.hbFrag2.max()),
        verdicts=sorted(set(d.E2_valence.map(verdict))), verdicts_withRY_vs_05=sorted(set((d.E2_with_RY >= 0.5).map({True: 'pass', False: 'fail'})))))
pt = pd.DataFrame(rows); pt.to_csv(os.path.join(V, 'protein_contacts_table.csv'), index=False)
S['protein'] = dict(contacts=pt.to_dict(orient='records'))
# QM size and environment effects (valence)
def span(d, vs): x = d[d.variant.isin(vs)]; return x.groupby('contact').E2_valence.agg(lambda s: s.max() - s.min()).max()
S['protein']['size_span'] = r2(span(pr, ['cluster_pcm4', 'shell3', 'shell4', 'shell5'])); S['protein']['ff_span'] = r2(span(pr, ['qmmm_ff14SB', 'qmmm_charmm36']))
S['protein']['boundary_span'] = r2(span(pr, ['qmmm_ff14SB', 'qmmm_Z1', 'qmmm_Z2', 'qmmm_Z3'])); S['protein']['radius_span'] = r2(span(pr, ['qmmm_r10', 'qmmm_ff14SB', 'qmmm_r20', 'qmmm_r25', 'qmmm_full']))
S['protein']['pcm_gas_span'] = r2(span(pr, ['cluster_pcm4', 'cluster_gas']))
# M12: change cluster -> QM/MM against the MM field E_AX
pm = pd.read_csv(os.path.join(R, 'protein_models.csv'))
m12 = []
for c, d in pr.groupby('contact'):
    if len(d) < 10: continue
    cl = d[d.variant == 'cluster_pcm4'].E2_valence.iloc[0]; cg = d[d.variant == 'cluster_gas'].E2_valence.iloc[0]
    for v in ['qmmm_ff14SB', 'qmmm_charmm36', 'qmmm_Z1', 'qmmm_Z2', 'qmmm_Z3', 'qmmm_r10', 'qmmm_r20', 'qmmm_r25', 'qmmm_full']:
        q = d[d.variant == v]
        if not len(q): continue
        f = pm[pm.model == q.name.iloc[0]].field_AX_VperA
        m12.append(dict(contact=c, variant=v, E_AX=float(f.iloc[0]) if len(f) else np.nan, dE_vs_gas=q.E2_valence.iloc[0] - cg, dE_vs_pcm=q.E2_valence.iloc[0] - cl, val=q.E2_valence.iloc[0]))
m12 = pd.DataFrame(m12)
# note: the field is computed at the first contact of each model; use only rows whose contact is that first contact
first = {r.model: r.xbs.split(' ; ')[0].replace(' ', '') for r in pm.itertuples() if isinstance(r.xbs, str)}
m12['first'] = [first.get(pr[(pr.contact == r.contact) & (pr.variant == r.variant)].name.iloc[0], '') == r.contact.replace(' ', '') for r in m12.itertuples()]
m12f = m12[m12['first'] & m12.E_AX.notna()]
d15 = m12f[m12f.variant == 'qmmm_ff14SB']
S['m12'] = dict(n_contacts=int(len(d15)), rho_across_default=rho(d15.E_AX, d15.dE_vs_gas), rows_default=d15[['contact', 'E_AX', 'dE_vs_gas', 'dE_vs_pcm']].round(3).values.tolist(),
    rho_all_variants=rho(m12f.E_AX, m12f.dE_vs_gas), n_all=int(len(m12f)),
    within=[(c, rho(g.E_AX, g.val), int(len(g))) for c, g in m12f.groupby('contact') if len(g) >= 5])
m12.to_csv(os.path.join(V, 'm12_field_vs_dE.csv'), index=False)
# ------------------------------------------------------------------ in-cluster pairs and amide lone pairs
pp = pd.read_csv(os.path.join(V, 'pairs_v3.csv')); am = pd.read_csv(os.path.join(V, 'amideN_v3.csv'))
S['pairs'] = dict(n=int(len(pp)), t1pass=int(pp.pass_T1.sum()), t1fail=int((~pp.pass_T1).sum()), fail_max_val=r2(pp[~pp.pass_T1].E2_valence.max()), fail_hb_max=r2(pp[~pp.pass_T1].hbY.max()),
    fail_rows=pp[~pp.pass_T1][['name', 'X', 'Y', 'd', 'theta1', 'E2_valence', 'hbY']].round(2).values.tolist(), thr113_hbY=r2(pp[pp.Y.str.contains('THR113')].hbY.max()))
S['amide'] = dict(n=int(len(am)), with_LP=int(am.lp_N.sum()), max_to_CX=r2(am.max_LP_N_to_sigmaCX.max()))
# ------------------------------------------------------------------ field (no QM)
fd = pd.read_csv(os.path.join(V, 'field_v3.csv')); fd['radius'] = fd.radius.astype(str)
def fval(**k):
    q = fd.copy()
    for a, b in k.items(): q = q[q[a].astype(str) == str(b)]
    return q.groupby('system').E_AX.first()
base15 = fval(ff='ff14SB', boundary='RC', radius=15, explicitH=0)
S['field'] = dict(explicitH_shift={s: r2(v) for s, v in (fval(ff='ff14SB', boundary='RC', radius=15, explicitH=6) - base15).items()},
    RCD_shift={s: r2(v) for s, v in (fval(ff='ff14SB', boundary='RCD', radius=15, explicitH=0) - base15).items()},
    r15_vs_full={s: r2(v) for s, v in (base15 - fval(ff='ff14SB', boundary='RC', radius='full', explicitH=0)).items()},
    r25_vs_full={s: r2(v) for s, v in (fval(ff='ff14SB', boundary='RC', radius=25, explicitH=0) - fval(ff='ff14SB', boundary='RC', radius='full', explicitH=0)).items()},
    r30_vs_full={s: r2(v) for s, v in (fval(ff='ff14SB', boundary='RC', radius=30, explicitH=0) - fval(ff='ff14SB', boundary='RC', radius='full', explicitH=0)).items()},
    r40_vs_full={s: r2(v) for s, v in (fval(ff='ff14SB', boundary='RC', radius=40, explicitH=0) - fval(ff='ff14SB', boundary='RC', radius='full', explicitH=0)).items()},
    nZero=fd[(fd.radius == '15') & (fd.boundary == 'RC') & (fd.explicitH == 0)].groupby(['system', 'ff'])[['nMM', 'nZero', 'nCharges', 'net']].first().reset_index().values.tolist(),
    full_default={s: r2(v) for s, v in fval(ff='ff14SB', boundary='RC', radius='full', explicitH=6).items()})
# ------------------------------------------------------------------ protein BXQI (G and valence E; cluster EDA not used)
pb = []
for r in pt.itertuples():
    if r.val_default is None or not np.isfinite(r.val_default): continue
    pb.append(dict(contact=r.contact, G=r.G, E=r.val_default, verdict=verdict(r.val_default), BXQI=round(bxqi(r.G, r.val_default, None), 1)))
pb = pd.DataFrame(pb); order = {'supported': 0, 'inconclusive': 1, 'not supported': 2}
pb = pb.sort_values(by=['verdict', 'BXQI'], key=lambda s: s.map(order) if s.name == 'verdict' else -s); pb.to_csv(os.path.join(V, 'bxqi_protein_table.csv'), index=False)
S['bxqi_protein'] = pb.values.tolist()
# ------------------------------------------------------------------ protein EDA and counterpoise in kcal/mol
pe = pd.read_csv(os.path.join(R, 'qm_protein_eda.csv'))
S['protein_eda'] = [[r.sys, r.lig, r2(r.frz / KCAL), r2(r.pol / KCAL), r2(r.ct / KCAL), r2(r.total / KCAL), r2(r.chi_CT), r2(r.g16_cp_kJ / KCAL), r2(r.orca_cp_kJ / KCAL)] for r in pe.itertuples()]
S['protein_eda_orca_maxdiff_kcal'] = r2(max(abs(r.total - r.orca_cp_kJ) for r in pe.itertuples()) / KCAL)
cpp = pd.read_csv(os.path.join(R, 'orca_cp_protein.csv')) if os.path.exists(os.path.join(R, 'orca_cp_protein.csv')) else None
S['protein_cp_cols'] = list(cpp.columns) if cpp is not None else None
# ------------------------------------------------------------------ hIYD
hn = pd.read_csv(os.path.join(R, 'hiyd', 'hiyd_nbo_reanalysis.csv')); hn['v3'] = hn.E2_bondx_pair.map(verdict)
S['hiyd_nbo'] = dict(crosstab=pd.crosstab(hn.lewis, hn.v3).to_dict(), supported_lewis=hn[(hn.lewis == 'strong') & (hn.v3 == 'supported')][['pathway', 'residue', 'acceptor_atom', 'E2_bondx_pair', 'hb_LPI_to_YH']].values.tolist(),
    supported_poor=hn[(hn.lewis == 'weak') & (hn.v3 == 'supported')][['pathway', 'residue', 'acceptor_atom', 'E2_bondx_pair']].values.tolist(),
    backbone_N_LP=int(hn[hn.acceptor_class.str.startswith('backbone N')].LP_on_acceptor.sum()), backbone_N_n=int(hn.acceptor_class.str.startswith('backbone N').sum()),
    fmn_hb=r2(hn[hn.residue == 'FMN'].hb_LPI_to_YH.max()), ge05=int((hn.E2_bondx_pair >= 0.5).sum()), ge05_lewis=int(((hn.E2_bondx_pair >= 0.5) & (hn.lewis == 'strong')).sum()))
t2 = pd.read_csv(os.path.join(V, 'hiyd_t2_v3.csv'))
S['hiyd_t2'] = dict(rows=int(len(t2)), windows_with_contact=int(t2.groupby(['pathway', 'xi']).ngroups), occ_max=r2(t2.occupancy.max()), n_ge_010=int((t2.occupancy >= 0.10).sum()),
    pocket_occ_max=r2(t2[t2.xi <= 8.8].occupancy.max()), tau_range=[r2(t2.tau_ps.min()), r2(t2.tau_ps.max())], tau_median=r2(t2.tau_ps.median()), tau_max=r2(t2.tau_max_ps.max()),
    top=t2.sort_values('D', ascending=False).head(10)[['pathway', 'xi', 'acceptor', 'atom', 'lewis', 'occupancy', 'tau_ps', 'D', 'G_best']].values.tolist(), dt_ps=sorted(t2.dt_ps.unique().tolist()))
# T5 on the trajectory: best persistent Lewis-base contacts with a published NBO pose of the same residue, atom and pathway
def amber(acc): return acc.split(':')[0]
hn['key'] = hn.pathway + '|' + hn.residue_amber.str.extract(r'([A-Z]+)(\d+)').apply(lambda r: r[0] + str(int(r[1]) + 1), axis=1) + '|' + hn.acceptor_atom
t2['key'] = t2.pathway + '|' + t2.acceptor.str.split(':').str[0] + '|' + t2.atom
best = t2.sort_values('D', ascending=False).drop_duplicates('key').merge(hn[['key', 'E2_bondx_pair', 'xi']].rename(columns={'xi': 'xi_pub'}), on='key')
best['BXQI_GDE'] = [round(bxqi(r.G_best, r.E2_bondx_pair, None, w=dict(G=0.25, D=0.20, E=0.30)) if False else 0, 1) for r in best.itertuples()]
def bx3(G, D, E):
    t = dict(G=G, D=D, E=1 - np.exp(-max(0, E) / 2.0)); w = dict(G=0.25, D=0.20, E=0.30)
    if any(v <= 0 for v in t.values()): return 0.0
    return float(100 * np.exp(sum(w[k] / sum(w.values()) * np.log(t[k]) for k in t)))
best['BXQI_GDE'] = [round(bx3(r.G_best, r.D, r.E2_bondx_pair), 1) for r in best.itertuples()]; best['verdict'] = best.E2_bondx_pair.map(verdict)
best[['key', 'xi', 'occupancy', 'tau_ps', 'D', 'G_best', 'xi_pub', 'E2_bondx_pair', 'verdict', 'BXQI_GDE']].to_csv(os.path.join(V, 'hiyd_t5_table.csv'), index=False)
S['hiyd_t5'] = best[['key', 'xi', 'occupancy', 'D', 'G_best', 'E2_bondx_pair', 'verdict', 'BXQI_GDE']].round(3).values.tolist()
json.dump(S, open(os.path.join(V, 'summary_v3.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps(S, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))[:200])
