#!/usr/bin/env python3
"""28_v3_round2_analysis.py
Analysis added after the second review:
  - Sigma E(2) per unit binding energy for each acceptor class of HAL59 (bias of one threshold)
  - which data points set the extremes of the grey-zone spreads, and the zone without the Thr113 contact
  - the smallest supported benchmark value
  - additional job sets: hIYD hard negatives (set G), 4DST with the intact amide (set J), and the published hIYD
    backbone-N poses with an acetyl cap (set K); the updated hIYD decision counts
  - the competing-hydrogen-bond flag on all protein and hIYD contacts
  - the residence factor of D at the saved-frame spacing of the hIYD windows
Usage: python3 28_v3_round2_analysis.py [ROOT]
Output: results/v3/summary_round2.json, results/v3/class_ratio.csv, results/v3/hiyd_decisions_updated.csv
"""
import os, sys, json
import numpy as np, pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
R = os.path.join(ROOT, 'results'); V3 = os.path.join(R, 'v3')
LO, HI, HBMIN = 0.3, 0.8, 0.5
verdict = lambda v: 'supported' if v >= HI else 'inconclusive' if v >= LO else 'not supported'
r2 = lambda x: None if x is None else round(float(x), 2)
S = {}
nbo = pd.read_csv(os.path.join(V3, 'nbo_v3.csv'))

# ---------------------------------------------------------------- HAL59 acceptor classes
q = pd.read_csv(os.path.join(R, 'qm_hal59.csv')).merge(nbo[nbo.set == 'A_HAL59'][['name', 'E2_valence', 'verdict']], on='name')
xb = q[~q.control].copy(); xb['ratio'] = xb.E2_valence / xb.ref_kcal
cls = xb.groupby('acc_el').ratio.agg(['count', 'median', 'min', 'max'])
cls = cls.rename(index={'pi': 'ring'}).sort_values('median', ascending=False)
cls.round(2).to_csv(os.path.join(V3, 'class_ratio.csv'))
S['class_ratio'] = {k: dict(n=int(v['count']), median=r2(v['median']), min=r2(v['min']), max=r2(v['max'])) for k, v in cls.iterrows()}
med = cls['median']
S['class_ratio_note'] = dict(N_over_ring=r2(med['N'] / med['ring']), O_over_N=r2(med['O'] / med['N']),
                             BE_at_08_N=r2(HI / med['N']), BE_at_08_O=r2(HI / med['O']), BE_at_08_ring=r2(HI / med['ring']))
sup = xb[xb.verdict == 'supported']
S['supported_min'] = [sup.loc[sup.E2_valence.idxmin(), 'name'], r2(sup.E2_valence.min())]
S['supported_max'] = [sup.loc[sup.E2_valence.idxmax(), 'name'], r2(sup.E2_valence.max())]
ring = xb[xb.acc_el == 'pi'][['name', 'ref_kcal', 'E2_valence', 'verdict']]
S['ring_dimers'] = ring.round(2).values.tolist()

# ---------------------------------------------------------------- grey zone: which points set the extremes
mb = pd.read_csv(os.path.join(R, 'qm_method_basis.csv')).merge(nbo[nbo.set == 'D_method_basis'][['name', 'E2_valence']], on='name')
r0 = mb[mb.level == 'wB97X-D__def2-TZVP'].set_index('base').E2_valence; mb['rel'] = mb.E2_valence / mb.base.map(r0)
x = mb[mb.level != 'wB97X-D__def2-TZVP']
S['level_extremes'] = dict(min=[x.loc[x.rel.idxmin(), 'base'], x.loc[x.rel.idxmin(), 'level'], r2(x.rel.min())],
                           max=[x.loc[x.rel.idxmax(), 'base'], x.loc[x.rel.idxmax(), 'level'], r2(x.rel.max())])
pr = pd.read_csv(os.path.join(R, 'qm_protein.csv')).merge(nbo[nbo.set == 'F_protein'][['name', 'label', 'E2_valence']], on=['name', 'label'])
rows = []
for c, d in pr.groupby('contact'):
    if len(d) < 5: continue
    v0 = d[d.variant == 'cluster_pcm4'].E2_valence.iloc[0]
    if v0 < 0.2: continue
    rr = d.E2_valence / v0
    rows.append(dict(contact=c, n=len(d), v0=v0, rmin=rr.min(), vmin=d.variant.iloc[rr.values.argmin()], rmax=rr.max(), vmax=d.variant.iloc[rr.values.argmax()]))
M = pd.DataFrame(rows)
def zone(mm):
    lo, hi = mm.rmin.min(), mm.rmax.max(); clo, chi = x.rel.min() * lo, x.rel.max() * hi
    return dict(model=[r2(lo), r2(hi)], combined=[r2(clo), r2(chi)], zone=[r2(0.5 / chi), r2(0.5 / clo)])
S['greyzone_all'] = zone(M); S['greyzone_without_Thr113'] = zone(M[~M.contact.str.contains('THR113')])
S['model_extremes'] = dict(min=M.loc[M.rmin.idxmin(), ['contact', 'vmin']].tolist() + [r2(M.rmin.min())], max=M.loc[M.rmax.idxmax(), ['contact', 'vmax']].tolist() + [r2(M.rmax.max())])
mt = M[~M.contact.str.contains('THR113')]
S['model_extremes_without_Thr113'] = dict(max=mt.loc[mt.rmax.idxmax(), ['contact', 'vmax']].tolist() + [r2(mt.rmax.max())])

# ---------------------------------------------------------------- competing H-bond flag on the protein models and contacts
pn = nbo[nbo.set == 'F_protein'].copy(); pn['hbMax'] = pn[['hbY', 'hbFrag2']].max(axis=1)
pn['hbFlag'] = (pn.hbMax >= HBMIN) & (pn.hbMax > pn.E2_valence)
pn['contact'] = pn.label
fl = pn[pn.hbFlag].groupby('label').name.apply(lambda s: sorted(set(v.split('__')[-1] for v in s))).to_dict()
S['protein_hbflag_models'] = {k: v for k, v in fl.items()}
S['protein_hbflag_counts'] = {k: [int(pn[(pn.label == k)].hbFlag.sum()), int((pn.label == k).sum())] for k in fl}

# ---------------------------------------------------------------- additional job sets
ex = pd.read_csv(os.path.join(V3, 'extra_nbo.csv'))
ex['flags'] = ex['flags'].fillna('')
G = ex[ex.set == 'G_hard_negatives'].copy(); G['kind'] = np.where(G.role.str.startswith('control'), 'Lewis base', 'poor Lewis base')
S['setG'] = {k: dict(n=int(len(d)), verdicts=d.verdict.value_counts().to_dict(), hbflag=int(d['flags'].str.contains('hbond').sum()),
                     supported=d[d.verdict == 'supported'][['role', 'E2_valence']].round(2).values.tolist()) for k, d in G.groupby('kind')}
S['setG_all_geometry_ok'] = bool((G.geometry_ok.astype(str).str.lower() == 'true').all())
J = ex[ex.set == 'J_4dst_amide']; S['setJ'] = dict(n=int(len(J)), max_valence=r2(J.E2_valence.max()), max_hb=r2(J[['hbY', 'hbFrag2']].max().max()), lp=J.lpOnY.astype(str).tolist())
K = ex[ex.set == 'K_hiyd_amide'].copy()
K['pose'] = K.role.str.replace('published pose ', '').str.replace('; cap H -> C(=O)CH3', '', regex=False)
S['setK'] = dict(n=int(len(K)), verdicts=K.verdict.value_counts().to_dict(), max_valence=[K.loc[K.E2_valence.idxmax(), 'pose'], r2(K.E2_valence.max())],
                 hbflag=int(K['flags'].str.contains('hbond').sum()), hb_range=[r2(K[K['flags'].str.contains('hbond')].hbFrag2.min()), r2(K[K['flags'].str.contains('hbond')].hbFrag2.max())],
                 amine_supported_before=K[K.E_pub_gate >= HI][['pose', 'E_pub_gate', 'E2_valence']].round(2).values.tolist(),
                 incomplete_job=K[K.terminated.astype(str).str.lower() != 'true'].name.tolist())
H = ex[ex.set == 'H_marginal_levels']; S['setH'] = dict(n=int(len(H)), finished=int((H.terminated.astype(str).str.lower() == 'true').sum()))
S['setG_ranges'] = dict(intra=[r2(G['intra_O_sigmaCI'].min()), r2(G['intra_O_sigmaCI'].max())])

# ---------------------------------------------------------------- hIYD: decisions with the acetyl-capped backbone N poses
h = pd.read_csv(os.path.join(R, 'hiyd', 'hiyd_nbo_reanalysis.csv'))
h['hbFlag'] = (h.hb_LPI_to_YH >= HBMIN) & (h.hb_LPI_to_YH > h.E2_bondx_pair)
h['source'] = 'published model'; h['E2_used'] = h.E2_bondx_pair; h['verdict_used'] = h.verdict_v3; h['hb_used'] = h.hb_LPI_to_YH
bb = h.acceptor_class == 'backbone N-H'
key = K.pose.str.extract(r'^(\S+) N \(pathway (\w), xi ([\d.]+)\)')
K['residue'] = key[0]; K['pathway'] = key[1]; K['xi'] = key[2].astype(float)
for i in h[bb].index:
    m = K[(K.residue == h.at[i, 'residue']) & (K.pathway == h.at[i, 'pathway']) & (np.isclose(K.xi, h.at[i, 'xi']))]
    assert len(m) == 1, (h.at[i, 'residue'], h.at[i, 'pathway'], h.at[i, 'xi'])
    m = m.iloc[0]
    h.at[i, 'E2_used'] = m.E2_valence; h.at[i, 'verdict_used'] = m.verdict; h.at[i, 'hb_used'] = m.hbFrag2; h.at[i, 'source'] = 'acetyl cap (set K)'
    h.at[i, 'hbFlag'] = ('hbond' in m['flags'])
h.to_csv(os.path.join(V3, 'hiyd_decisions_updated.csv'), index=False)
S['hiyd_published_counts'] = h.verdict_v3.value_counts().to_dict()
S['hiyd_updated_counts'] = h.verdict_used.value_counts().to_dict()
S['hiyd_updated_supported'] = h[h.verdict_used == 'supported'][['pathway', 'residue', 'acceptor_atom', 'E2_used', 'hb_used', 'hbFlag']].round(2).values.tolist()
S['hiyd_updated_hbflag'] = h[h.hbFlag][['pathway', 'residue', 'acceptor_atom', 'E2_used', 'hb_used', 'verdict_used']].round(2).values.tolist()
S['hiyd_backbone_updated'] = h[bb].verdict_used.value_counts().to_dict()
S['hiyd_backbone_before'] = h[bb].verdict_v3.value_counts().to_dict()

# ---------------------------------------------------------------- residence factor at the hIYD frame spacing
t2 = pd.read_csv(os.path.join(V3, 'hiyd_t2_v3.csv'))
dts = sorted(set(t2.dt_ps.astype(int)))
S['hiyd_dt_with_contacts'] = dts
fac = lambda dt: (1 - np.exp(-dt / 100.0)) ** 0.25
S['residence_factor'] = {str(dt): r2(fac(dt)) for dt in [100] + dts}
S['hiyd_tau'] = dict(single_frame_share=r2((t2.tau_ps <= t2.dt_ps + 1e-6).mean()), n=int(len(t2)), median=r2(t2.tau_ps.median()), range=[r2(t2.tau_ps.min()), r2(t2.tau_ps.max())], max_run=r2(t2.tau_max_ps.max()), D_over_f_min=r2((t2.D / t2.occupancy).min()))

json.dump(S, open(os.path.join(V3, 'summary_round2.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps(S, indent=1, default=str)[:6000])
