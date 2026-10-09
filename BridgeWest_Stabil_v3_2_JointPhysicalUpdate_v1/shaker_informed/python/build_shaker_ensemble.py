from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)

fits = pd.read_csv(DATA / 'shaker_modal_fits.csv')
accepted = fits[fits['accepted'].astype(str).str.lower().eq('true')].copy()
d24 = accepted[accepted.year.eq(2024)].copy()
d25 = accepted[accepted.year.eq(2025)].copy()

# Nominal low-force state uses all accepted 5 m/s^2 calibration evidence.
# This is a central low-amplitude state, not a claim that the bridge has one
# campaign-invariant modal state.
base = {}
for mode in (1, 2):
    low = accepted[(accepted['mode'] == mode) & np.isclose(accepted.shaker_acc_rms_m_s2, 5.0)]
    base[mode] = dict(f_Hz=float(low.frequency_Hz.median()),
                      zeta=float(low.damping_ratio.median()), n=int(len(low)))

# 2024: genuine same-run joint states.
piv = d24.pivot_table(index=['run', 'shaker_acc_rms_m_s2'], columns='mode',
    values=['frequency_Hz','damping_ratio','MAC_attached','fit_relative_L2','window_sensitivity'], aggfunc='first')
rows24=[]
for (run,level),r in piv.iterrows():
    needed=[('frequency_Hz',1),('frequency_Hz',2),('damping_ratio',1),('damping_ratio',2)]
    if not all(k in r.index and np.isfinite(r[k]) for k in needed): continue
    f1,f2=float(r[('frequency_Hz',1)]),float(r[('frequency_Hz',2)])
    z1,z2=float(r[('damping_ratio',1)]),float(r[('damping_ratio',2)])
    rows24.append(dict(source='2024_same_run',year=2024,level_m_s2=float(level),
        source_mode1=f'run{int(run):02d}',source_mode2=f'run{int(run):02d}',
        f1_Hz=f1,f2_Hz=f2,zeta1=z1,zeta2=z2,
        MAC1=float(r.get(('MAC_attached',1),np.nan)),MAC2=float(r.get(('MAC_attached',2),np.nan))))
S24=pd.DataFrame(rows24)
S24.insert(0,'state_2024',np.arange(1,len(S24)+1))
S24.to_csv(OUT/'shaker_joint_states_2024.csv',index=False)

# 2025: modes were acquired in separate run blocks. Do not rank-pair them.
# For the production uncertainty mixture, represent their joint state as the
# product of the two empirical marginals *within the same excitation level*.
# This encodes independence conditional on the known shaker level rather than
# inventing same-run correlation. Campaigns receive equal total weight; within
# each campaign the three shaker levels receive equal total weight.
rows25=[]
for level in sorted(d25.shaker_acc_rms_m_s2.unique()):
    m1=d25[(d25['mode']==1)&np.isclose(d25.shaker_acc_rms_m_s2,level)]
    m2=d25[(d25['mode']==2)&np.isclose(d25.shaker_acc_rms_m_s2,level)]
    if m1.empty or m2.empty: continue
    nprod=len(m1)*len(m2)
    for _,r1 in m1.iterrows():
        for _,r2 in m2.iterrows():
            rows25.append(dict(source='2025_independent_within_level',year=2025,level_m_s2=float(level),
                source_mode1=f'run{int(r1.run):02d}',source_mode2=f'run{int(r2.run):02d}',
                f1_Hz=float(r1.frequency_Hz),f2_Hz=float(r2.frequency_Hz),
                zeta1=float(r1.damping_ratio),zeta2=float(r2.damping_ratio),
                MAC1=float(r1.MAC_attached),MAC2=float(r2.MAC_attached),
                component_weight=0.5/3.0/nprod))
S25=pd.DataFrame(rows25)
S25.insert(0,'state_2025_product',np.arange(1,len(S25)+1))
S25.to_csv(OUT/'shaker_product_states_2025.csv',index=False)

# Assign 2024 component weights after knowing paired counts at each level.
parts=[]
for level,d in S24.groupby('level_m_s2',sort=True):
    x=d.copy(); x['component_weight']=0.5/3.0/len(x); parts.append(x)
W24=pd.concat(parts,ignore_index=True)
W24=W24.drop(columns=['state_2024'])
W25=S25.drop(columns=['state_2025_product'])
W=pd.concat([W24,W25],ignore_index=True,sort=False)
W.insert(0,'state',np.arange(1,len(W)+1))
W['weight']=W['component_weight']
W=W.drop(columns=['component_weight'])
for mode in (1,2):
    W[f'lambda_scale_mode{mode}']=(W[f'f{mode}_Hz']/base[mode]['f_Hz'])**2
    W[f'zeta_scale_mode{mode}']=W[f'zeta{mode}']/base[mode]['zeta']
assert abs(W.weight.sum()-1)<1e-12
W.to_csv(OUT/'shaker_weighted_states_all.csv',index=False)

# 2025 single-fit cross-campaign audit against same-level 2024 marginals.
check=[]
for _,r in d25.iterrows():
    mode=int(r['mode']); level=float(r.shaker_acc_rms_m_s2)
    dev=d24[(d24['mode']==mode)&np.isclose(d24.shaker_acc_rms_m_s2,level)]
    fm,fs=float(dev.frequency_Hz.mean()),float(dev.frequency_Hz.std(ddof=1))
    zm,zs=float(dev.damping_ratio.mean()),float(dev.damping_ratio.std(ddof=1))
    check.append(dict(run=int(r.run),year=2025,mode=mode,shaker_acc_rms_m_s2=level,
        frequency_Hz=float(r.frequency_Hz),damping_ratio=float(r.damping_ratio),
        MAC_attached=float(r.MAC_attached),fit_relative_L2=float(r.fit_relative_L2),
        window_sensitivity=float(r.window_sensitivity),
        z_frequency_vs_2024=(float(r.frequency_Hz)-fm)/fs if fs>0 else np.nan,
        z_damping_vs_2024=(float(r.damping_ratio)-zm)/zs if zs>0 else np.nan))
pd.DataFrame(check).to_csv(OUT/'shaker_crosscampaign_2025_unpaired.csv',index=False)

# Weighted prior moments of the modal ensemble.
x=W[['f1_Hz','f2_Hz','zeta1','zeta2']].to_numpy(); w=W.weight.to_numpy(); mu=(w[:,None]*x).sum(0)
xc=x-mu; cov=xc.T@(xc*w[:,None]); sd=np.sqrt(np.diag(cov)); corr=cov/np.outer(sd,sd)
summary={
 'archive_records_total':113,
 'accepted_single_mode_fits_total':int(len(accepted)),
 'accepted_single_mode_fits_2024':int(len(d24)),
 'accepted_single_mode_fits_2025':int(len(d25)),
 'genuine_joint_states_2024':int(len(S24)),
 'independent_product_states_2025':int(len(S25)),
 'production_weighted_states':int(len(W)),
 'production_weighting':'0.5 total weight per shaker campaign; within each campaign 1/3 per excitation level; 2024 uses genuine same-run pairs; 2025 uses the product of empirical mode marginals within level because modes were measured in separate blocks.',
 'baseline_policy':'Median of all accepted 5 m/s^2 fits for each mode. The weighted ensemble represents amplitude and campaign variation around this central low-force state.',
 'baseline':{'mode1_frequency_Hz':base[1]['f_Hz'],'mode2_frequency_Hz':base[2]['f_Hz'],
             'mode1_damping_ratio':base[1]['zeta'],'mode2_damping_ratio':base[2]['zeta'],
             'n_mode1_lowforce':base[1]['n'],'n_mode2_lowforce':base[2]['n']},
 'weighted_mean':dict(zip(['f1_Hz','f2_Hz','zeta1','zeta2'],mu.tolist())),
 'weighted_covariance':cov.tolist(),'weighted_correlation':corr.tolist(),
 'ranges_all_weighted_states':{k:[float(W[k].min()),float(W[k].max())] for k in ['f1_Hz','f2_Hz','zeta1','zeta2','lambda_scale_mode1','lambda_scale_mode2']},
 'stabil_mode_mapping':{'first_global_bending_fe_index':3,'first_global_torsion_fe_index':8,
   'nominal_rebuild_frequency_Hz':[5.211011589665426,7.736400263217337],
   'MAC_to_2025_four_station_shaker_shapes':[0.998747,0.996079],
   'policy':'Mode shapes are retained; the two identified global modal eigenvalues and Rayleigh damping anchors are updated to each shaker state.'},
 'interpretation':'A weighted empirical/modal uncertainty measure for discrepancy propagation. It is not a probability distribution of Young modulus or unique physical section stiffness.'
}
(OUT/'shaker_ensemble_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
