from __future__ import annotations
from pathlib import Path
import json, sys, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.optimize import minimize
from scipy.special import ndtr, ndtri, log_ndtr, ndtri_exp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]; SROOT=ROOT/'shaker_informed'; OUT=SROOT/'results'/'run40_continuous_load'
OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(HERE))
from stabil_shell_observation import build_bridgewest_strain_modal, self_test
from run40_shaker_informed_virtual_strain import modal_parameters, force_history, simulate_basis, weighted_discrepancy, scores, CHANNELS, REF, TARGET, TARGET_NAMES

PROCESSED=Path('/mnt/data/bridgewest_work/restored/BridgeWest_traffic_simulation_restored/inputs/processed/2025_run40.mat')


def _rtruncstd(a,b,rng):
    # Stable scalar N(0,1) draw on [a,b]. Positive-tail intervals are
    # mirrored so the computation is performed in the negative tail.
    if a>0:
        return -_rtruncstd(-b,-a,rng)
    if b>0:
        Fa=ndtr(a);Fb=ndtr(b)
        u=Fa+(Fb-Fa)*rng.random()
        return ndtri(np.clip(u,1e-300,1-1e-16))
    la=log_ndtr(a);lb=log_ndtr(b)
    r=np.exp(la-lb) if np.isfinite(la) else 0.0
    q=r+(1-r)*rng.random()
    return ndtri_exp(lb+np.log(q))

def bounded_map(H,h,lo,hi):
    mid=(lo+hi)/2
    def f(x): return .5*x@H@x-h@x
    def g(x): return H@x-h
    r=minimize(f,mid,jac=g,bounds=list(zip(lo,hi)),method='L-BFGS-B',options={'ftol':1e-14,'gtol':1e-10,'maxiter':2000})
    if not r.success: print('MAP warning:',r.message)
    return r.x


def truncated_gaussian_moments(H,h,lo,hi,seed=20260924,burn=1500,keep=8000,thin=1):
    H=(H+H.T)/2
    ev=np.linalg.eigvalsh(H); ridge=max(0.0,1e-12*max(ev.max(),1.0)-ev.min())
    if ridge>0:H=H+np.eye(len(lo))*ridge
    m=np.linalg.solve(H,h);x=bounded_map(H,h,lo,hi);rng=np.random.default_rng(seed)
    out=np.empty((keep,len(lo)));n=0;total=burn+keep*thin
    for sweep in range(total):
        for i in range(len(lo)):
            hii=H[i,i];sd=np.sqrt(1/hii)
            cm=m[i]-(H[i]@(x-m)-hii*(x[i]-m[i]))/hii
            a=(lo[i]-cm)/sd;b=(hi[i]-cm)/sd
            x[i]=cm+sd*_rtruncstd(a,b,rng)
        if sweep>=burn and (sweep-burn)%thin==0:
            out[n]=x;n+=1
    mu=out.mean(0);C=np.cov(out,rowvar=False,ddof=1)
    return mu,C,out,m,ev


def load_posterior(Bref_kN,yref,mu_ref,invR,ti,lo,hi,seed):
    H=np.zeros((8,8));h=np.zeros(8)
    for t in ti:
        B=Bref_kN[t];S=invR[t];d=yref[t]-mu_ref[t]
        H+=B.T@S@B;h+=B.T@S@d
    mu,C,samp,uncon,ev=truncated_gaussian_moments(H,h,lo,hi,seed=seed)
    return mu,C,samp,H,h,uncon,ev


def main():
    start=time.time();assert self_test()<1e-12
    fe=loadmat(ROOT/'results'/'python_model.mat',simplify_cells=True);cfg=json.loads(str(fe['config_json']))
    nodes=np.asarray(fe['nodes'],float);shells=np.asarray(fe['shells'],float);P=np.asarray(fe['modes'],float);freq0=np.asarray(fe['frequencies_Hz'],float).ravel()
    summ=json.loads((SROOT/'results'/'shaker_ensemble_summary.json').read_text());states=pd.read_csv(SROOT/'results'/'shaker_weighted_states_all.csv');base=summ['baseline']
    ib=summ['stabil_mode_mapping']['first_global_bending_fe_index']-1;it=summ['stabil_mode_mapping']['first_global_torsion_fe_index']-1
    sel=np.where(freq0<=25)[0];sel=np.unique(np.r_[sel,ib,it]);Psel=P[:,sel]
    Bstrain,audit=build_bridgewest_strain_modal(nodes,shells,P,cfg,sel);pd.DataFrame(audit).to_csv(OUT/'strain_observation_map_audit.csv',index=False)
    z=loadmat(PROCESSED,simplify_cells=True);tobs=np.asarray(z['time'],float).ravel();y=np.asarray(z['Y'],float);speed=float(z['speed'])
    dt=.0025;ti_internal=np.arange(tobs[0],tobs[-1]+dt/2,dt);sample_idx=np.rint((tobs-tobs[0])/dt).astype(int)
    ax=np.array([0,2.6,19,21.6,25.9,28.5,44.9,47.5,51.8,54.4,70.8,73.4,77.7,80.3,96.7,99.3])
    pm=force_history(nodes,shells,Psel,sel,ti_internal,speed,np.asarray(cfg['bearing_x_m'],float),ax,float(np.mean(cfg['girder_y_m'])))
    om,zeta,ob,zb=modal_parameters(freq0,sel,ib,it,states,base)
    B0=simulate_basis(pm,Bstrain,ob[None,:],zb[None,:],dt,sample_idx)[0].astype(float)*1e6 # microstrain/N
    Bs=simulate_basis(pm,Bstrain,om,zeta,dt,sample_idx).astype(float)*1e6
    DB=Bs-B0[None,:,:,:]
    sw=states.weight.to_numpy();sw=sw/sw.sum()
    c=dict(cov_ridge=.05,ref_noise=.5,target_noise=.5,likelihood_stride=10,interval_z=1.959964)
    lidx=np.arange(0,len(tobs),c['likelihood_stride'])
    lo=np.full(8,180.0);hi=np.full(8,360.0);pbar=np.full(8,270.0)
    trace=[]
    for itn in range(20):
        delta=np.einsum('stob,b->sto',DB,pbar*1e3,optimize=True)
        mu,gain,cv,invR=weighted_discrepancy(delta,sw,REF,c)
        Bref=B0[:,REF,:]*1e3 # µe/kN
        pmean,Pcov,samples,H,h,uncon,ev=load_posterior(Bref,y[:,REF],mu[:,REF],invR,lidx,lo,hi,seed=20260924+itn)
        change=np.max(np.abs(pmean-pbar));trace.append(dict(iteration=itn+1,max_bogie_change_kN=float(change),condition_number=float(np.linalg.cond(H)),minimum_precision_eigenvalue=float(ev.min()),maximum_precision_eigenvalue=float(ev.max()),**{f'p{i+1}_kN':float(pmean[i]) for i in range(8)}))
        pbar=pmean
        if change<.05:break
    pd.DataFrame(trace).to_csv(OUT/'reference_only_iteration_trace.csv',index=False)
    # Analytic affine propagation of truncated load posterior moments.
    T,O,_=B0.shape;mean=np.zeros((T,O));sd=np.zeros((T,O))
    for t in range(T):
        B=B0[t]*1e3;Br=B[REF];G=gain[t]
        A=B-G@Br
        cc=mu[t]+G@(y[t,REF]-mu[t,REF])
        mean[t]=cc+A@pmean
        sd[t]=np.sqrt(np.maximum(np.einsum('oi,ij,oj->o',A,Pcov,A,optimize=True)+cv[t],0))
    low=mean-c['interval_z']*sd;high=mean+c['interval_z']*sd
    sc=scores(y,mean,low,high);sc.to_csv(OUT/'run40_target_scores.csv',index=False)
    # Save posterior load moments/correlations.
    ltab=pd.DataFrame({'bogie':np.arange(1,9),'posterior_mean_kN':pmean,'posterior_sd_kN':np.sqrt(np.diag(Pcov)),'lower_bound_kN':lo[:8] if False else np.full(8,180.),'upper_bound_kN':np.full(8,360.)})
    ltab.to_csv(OUT/'run40_posterior_loads.csv',index=False)
    np.savetxt(OUT/'run40_load_posterior_correlation.csv',Pcov/np.sqrt(np.outer(np.diag(Pcov),np.diag(Pcov))),delimiter=',')
    # protected histories opened only for assessment here
    hist={'time_s':tobs}
    for i,name in enumerate(CHANNELS):
        hist[f'{name}_measured']=y[:,i];hist[f'{name}_mean']=mean[:,i];hist[f'{name}_lo']=low[:,i];hist[f'{name}_hi']=high[:,i]
    pd.DataFrame(hist).to_csv(OUT/'run40_histories.csv',index=False)
    refrows=[]
    for idx in REF:
        rng=np.ptp(y[:,idx]);refrows.append(dict(channel=CHANNELS[idx],nrmse=float(np.sqrt(np.mean((mean[:,idx]-y[:,idx])**2))/rng),correlation=float(np.corrcoef(mean[:,idx],y[:,idx])[0,1]),coverage=float(np.mean((y[:,idx]>=low[:,idx])&(y[:,idx]<=high[:,idx])))))
    pd.DataFrame(refrows).to_csv(OUT/'run40_reference_fit.csv',index=False)
    # plots
    fig,axs=plt.subplots(2,2,figsize=(11,7),constrained_layout=True);axs=axs.ravel()
    for ap,idx,name in zip(axs,TARGET,TARGET_NAMES):
        ap.fill_between(tobs,low[:,idx],high[:,idx],alpha=.2,label='Nominal 95% interval')
        ap.plot(tobs,y[:,idx],label='Measured withheld');ap.plot(tobs,mean[:,idx],label='Prediction')
        row=sc.loc[sc.target==name].iloc[0];ap.set_title(f'{name}: NRMSE={row.nrmse:.3f}, r={row.correlation:.3f}')
        ap.set_xlabel('Time from bridge entry (s)');ap.set_ylabel('Strain (µε)');ap.grid(True,alpha=.2)
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,bbox_to_anchor=(.5,-.02));fig.savefig(OUT/'run40_targets.png',dpi=220,bbox_inches='tight');plt.close(fig)
    fig,axp=plt.subplots(figsize=(8.5,4.5),constrained_layout=True);x=np.arange(1,9);axp.errorbar(x,pmean,yerr=np.sqrt(np.diag(Pcov)),fmt='o',capsize=3);axp.axhline(180,linestyle='--');axp.axhline(360,linestyle='--');axp.set_xticks(x);axp.set_xlabel('Bogie');axp.set_ylabel('Load (kN)');axp.grid(True,alpha=.2);fig.savefig(OUT/'run40_load_posterior.png',dpi=220);plt.close(fig)
    out=dict(event='2025_run40',method='bounded continuous Gaussian posterior for 8 bogie loads + 84-state shaker-informed Stabil discrepancy',speed_m_s=speed,shaker_states=len(states),shaker_effective_state_count=float(1/np.sum(sw**2)),selected_dynamic_modes=len(sel),likelihood_samples=len(lidx),load_bounds_kN=[180,360],posterior_mean_bogie_loads_kN=pmean.tolist(),posterior_sd_bogie_loads_kN=np.sqrt(np.diag(Pcov)).tolist(),load_precision_condition_number=float(np.linalg.cond(H)),iteration_count=len(trace),final_max_bogie_change_kN=trace[-1]['max_bogie_change_kN'],targets=sc.to_dict(orient='records'),references=refrows,elapsed_s=time.time()-start,caveat='Run 40 was previously inspected and is a development diagnostic, not sealed prospective validation. Withheld target histories enter only the final scoring/plotting stage.')
    (OUT/'run40_summary.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

if __name__=='__main__':main()
