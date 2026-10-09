from pathlib import Path
import json,sys
import numpy as np,pandas as pd
from scipy.io import loadmat
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];SROOT=ROOT/'shaker_informed';OUT=SROOT/'results'/'run40_continuous_load'
sys.path.insert(0,str(HERE))
from stabil_shell_observation import build_bridgewest_strain_modal,self_test
from run40_shaker_informed_virtual_strain import modal_parameters,force_history,simulate_basis,CHANNELS,REF,TARGET,TARGET_NAMES
PROCESSED=Path('/mnt/data/bridgewest_work/restored/BridgeWest_traffic_simulation_restored/inputs/processed/2025_run40.mat')

def main():
 assert self_test()<1e-12
 fe=loadmat(ROOT/'results'/'python_model.mat',simplify_cells=True);cfg=json.loads(str(fe['config_json']))
 nodes=np.asarray(fe['nodes'],float);shells=np.asarray(fe['shells'],float);P=np.asarray(fe['modes'],float);freq0=np.asarray(fe['frequencies_Hz'],float).ravel()
 summ=json.loads((SROOT/'results'/'shaker_ensemble_summary.json').read_text());states=pd.read_csv(SROOT/'results'/'shaker_weighted_states_all.csv');base=summ['baseline']
 ib=summ['stabil_mode_mapping']['first_global_bending_fe_index']-1;it=summ['stabil_mode_mapping']['first_global_torsion_fe_index']-1
 sel=np.unique(np.r_[np.where(freq0<=25)[0],ib,it]);Psel=P[:,sel];Bstrain,_=build_bridgewest_strain_modal(nodes,shells,P,cfg,sel)
 z=loadmat(PROCESSED,simplify_cells=True);tobs=np.asarray(z['time'],float).ravel();y=np.asarray(z['Y'],float);speed=float(z['speed'])
 dt=.0025;ti=np.arange(tobs[0],tobs[-1]+dt/2,dt);sample_idx=np.rint((tobs-tobs[0])/dt).astype(int)
 ax=np.array([0,2.6,19,21.6,25.9,28.5,44.9,47.5,51.8,54.4,70.8,73.4,77.7,80.3,96.7,99.3])
 pm=force_history(nodes,shells,Psel,sel,ti,speed,np.asarray(cfg['bearing_x_m'],float),ax,float(np.mean(cfg['girder_y_m'])))
 om,zeta,ob,zb=modal_parameters(freq0,sel,ib,it,states,base)
 B0=simulate_basis(pm,Bstrain,ob[None,:],zb[None,:],dt,sample_idx)[0].astype(float)*1e6
 Bs=simulate_basis(pm,Bstrain,om,zeta,dt,sample_idx).astype(float)*1e6;DB=Bs-B0[None]
 loads=pd.read_csv(OUT/'run40_posterior_loads.csv').posterior_mean_kN.to_numpy()*1e3
 d=np.einsum('stob,b->sto',DB,loads,optimize=True);w=states.weight.to_numpy();w=w/w.sum();mu=np.einsum('s,sto->to',w,d)
 dc=d-mu[None];var=np.einsum('s,sto,sto->to',w,dc,dc,optimize=True)
 # conditional variance and state correlations, without measurement noise in reduction fraction
 rows=[]
 nom=np.einsum('tob,b->to',B0,loads,optimize=True)
 active=np.max(np.abs(nom),axis=1)>.05*np.max(np.abs(nom))
 for idx,name in zip(TARGET,TARGET_NAMES):
  cond=[];corrs=[[],[]]
  for t in range(len(tobs)):
   D=dc[:,t,REF];C=D.T@(D*w[:,None]);cross=dc[:,t,idx]@(D*w[:,None]);
   Sinv=np.linalg.pinv(C+np.eye(2)*(.05*max(np.trace(C)/2,1e-12)+1e-12+.5**2),rcond=1e-15)
   cv=max(var[t,idx]-cross@Sinv@cross,0)
   cond.append(cv)
   for jj,r in enumerate(REF):
    den=np.sqrt(max(var[t,idx]*var[t,r],1e-30));corrs[jj].append((dc[:,t,idx]@(dc[:,t,r]*w))/den)
  cond=np.array(cond); corr=np.array(corrs).T
  q=active & (var[:,idx]>1e-10)
  residual=y[:,idx]-(nom[:,idx]+mu[:,idx])
  rows.append(dict(target=name,
    rms_unconditioned_shaker_sd_microstrain=float(np.sqrt(np.mean(var[q,idx]))),
    rms_conditional_shaker_sd_microstrain=float(np.sqrt(np.mean(cond[q]))),
    median_variance_reduction_fraction=float(np.median(1-cond[q]/np.maximum(var[q,idx],1e-30))),
    median_abs_corr_with_A031=float(np.median(np.abs(corr[q,0]))),
    median_abs_corr_with_B053=float(np.median(np.abs(corr[q,1]))),
    rms_nominal_plus_shakermean_residual_microstrain=float(np.sqrt(np.mean(residual[q]**2)))))
 df=pd.DataFrame(rows);df.to_csv(OUT/'run40_shaker_observability.csv',index=False)
 fig,axp=plt.subplots(figsize=(8.5,4.8),constrained_layout=True);x=np.arange(len(df));width=.36
 axp.bar(x-width/2,df.rms_unconditioned_shaker_sd_microstrain,width,label='Unconditioned shaker SD')
 axp.bar(x+width/2,df.rms_conditional_shaker_sd_microstrain,width,label='After A031/B053 conditioning')
 axp.set_xticks(x,df.target);axp.set_ylabel('RMS shaker discrepancy SD (µε)');axp.grid(True,axis='y',alpha=.2);axp.legend()
 fig.savefig(OUT/'run40_shaker_observability.png',dpi=220);plt.close(fig)
 print(df.to_string(index=False))
if __name__=='__main__':main()
