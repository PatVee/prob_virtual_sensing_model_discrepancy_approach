from __future__ import annotations
from pathlib import Path
import json, sys, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SROOT=ROOT/'shaker_informed'
sys.path.insert(0,str(HERE))
from stabil_shell_observation import build_bridgewest_strain_modal, self_test

# Mounted current-session inputs from the earlier Bridge West traffic package.
TRAFFIC=Path('/mnt/data/bridgewest_work/restored/BridgeWest_traffic_simulation_restored/BridgeWest_v2_MATLAB_four_targets')
PROCESSED=Path('/mnt/data/bridgewest_work/restored/BridgeWest_traffic_simulation_restored/inputs/processed/2025_run40.mat')
DRAWS=TRAFFIC/'data'/'2025_run40_draws.mat'
OUT=SROOT/'results'/'run40_shaker_informed'
OUT.mkdir(parents=True,exist_ok=True)

CHANNELS=['A031','B031','A053_lower','B053','A058','A070','B081']
REF=np.array([0,3]); TARGET=np.array([4,5,1,6])
TARGET_NAMES=[CHANNELS[i] for i in TARGET]


def modal_parameters(freq0,sel,ib,it,states,base):
    nstate=len(states);freq=np.tile(freq0[sel],(nstate,1))
    posb=int(np.flatnonzero(sel==ib)[0]);post=int(np.flatnonzero(sel==it)[0])
    freq[:,posb]=states.f1_Hz.to_numpy();freq[:,post]=states.f2_Hz.to_numpy()
    om=2*np.pi*freq
    z=np.empty_like(om)
    for s,r in enumerate(states.itertuples(index=False)):
        o1=2*np.pi*r.f1_Hz;o2=2*np.pi*r.f2_Hz
        A=np.array([[1/(2*o1),o1/2],[1/(2*o2),o2/2]])
        alpha,beta=np.linalg.solve(A,np.array([r.zeta1,r.zeta2]))
        z[s]=np.clip(alpha/(2*om[s])+beta*om[s]/2,.002,.08)
        z[s,posb]=r.zeta1;z[s,post]=r.zeta2
    # nominal state
    fbase=freq0[sel].copy();fbase[posb]=base['mode1_frequency_Hz'];fbase[post]=base['mode2_frequency_Hz']
    ob=2*np.pi*fbase;o1=2*np.pi*base['mode1_frequency_Hz'];o2=2*np.pi*base['mode2_frequency_Hz']
    A=np.array([[1/(2*o1),o1/2],[1/(2*o2),o2/2]])
    alpha,beta=np.linalg.solve(A,np.array([base['mode1_damping_ratio'],base['mode2_damping_ratio']]))
    zb=np.clip(alpha/(2*ob)+beta*ob/2,.002,.08);zb[posb]=base['mode1_damping_ratio'];zb[post]=base['mode2_damping_ratio']
    return om,z,ob,zb


def force_history(nodes,shells,Psel,sel,time_internal,speed,bearing,axles,ytrack):
    phi_z=Psel[2::6,:] if Psel.shape[0] == 6*len(nodes) else None
    # Psel passed as full dof x selected modes.
    phi_z=Psel[2::6,:]
    deck=np.flatnonzero(shells[:,5].astype(int)==4)
    # Track-centreline deck cells, sorted by x minimum.
    cells=[]
    for e in deck:
        ids=shells[e,:4].astype(int);p=nodes[ids]
        if ytrack>=p[:,1].min()-1e-10 and ytrack<=p[:,1].max()+1e-10:
            cells.append((p[:,0].min(),p[:,0].max(),e))
    cells.sort(key=lambda x:(x[0],x[1]))
    xmins=np.array([c[0] for c in cells]);xmaxs=np.array([c[1] for c in cells])
    pm=np.zeros((len(time_internal),len(sel),8),dtype=float)
    a,b=bearing
    last_e=None;last_data=None
    for it,t in enumerate(time_internal):
        xlead=a+speed*t
        for ia,off in enumerate(axles):
            xp=xlead-off
            if xp<a-1e-12 or xp>b+1e-12:continue
            # choose first centreline cell containing x; boundary duplicates are harmless.
            candidates=np.flatnonzero((xmins<=xp+1e-10)&(xmaxs>=xp-1e-10))
            if len(candidates)==0:continue
            e=cells[candidates[0]][2]
            if e!=last_e:
                ids=shells[e,:4].astype(int);p=nodes[ids]
                last_data=(ids,p,p[:,0].min(),p[:,0].max(),p[:,1].min(),p[:,1].max());last_e=e
            ids,p,xmin,xmax,ymin,ymax=last_data
            xi=2*(xp-xmin)/(xmax-xmin)-1;eta=2*(ytrack-ymin)/(ymax-ymin)-1
            N=.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)])
            pm[it,:,ia//2] += -0.5*(N[:,None]*phi_z[ids]).sum(0)
    return pm


def simulate_basis(pm,obs_modal,om,zeta,dt,sample_idx):
    """Vectorized Newmark over multiple structural states.
    om/zeta are state x mode. Returns state x Tobs x output x 8.
    """
    ns,nm=om.shape;nb=pm.shape[2];no=obs_modal.shape[0]
    kk=om**2;cc=2*zeta*om
    beta=.25;gamma=.5;a0=1/(beta*dt**2);a1=gamma/(beta*dt);a2=1/(beta*dt);a3=1/(2*beta)-1;a4=gamma/beta-1;a5=dt*(gamma/(2*beta)-1)
    den=kk+a0+a1*cc
    q=np.zeros((ns,nm,nb));qd=np.zeros_like(q);qdd=np.zeros_like(q)
    out=np.empty((ns,len(sample_idx),no,nb),dtype=np.float32);isamp=0
    sample_set={int(k):j for j,k in enumerate(sample_idx)}
    for j in range(pm.shape[0]):
        force=pm[j][None,:,:]
        if j==0:
            qdd=force-cc[:,:,None]*qd-kk[:,:,None]*q
        else:
            pe=force+a0*q+a2*qd+a3*qdd+cc[:,:,None]*(a1*q+a4*qd+a5*qdd)
            qn=pe/den[:,:,None];an=a0*(qn-q)-a2*qd-a3*qdd;vn=qd+dt*((1-gamma)*qdd+gamma*an)
            q,qd,qdd=qn,vn,an
        if j in sample_set:
            jj=sample_set[j]
            # state x output x bogie
            out[:,jj]=np.einsum('om,smb->sob',obs_modal,q,optimize=True)
    return out


def weighted_discrepancy(delta,state_w,refs,c):
    # delta state x time x output
    sw=state_w/state_w.sum();mu=np.einsum('s,sto->to',sw,delta,optimize=True)
    ns,T,O=delta.shape;R=len(refs)
    gain=np.zeros((T,O,R));cv=np.zeros((T,O));invR=np.zeros((T,R,R))
    for t in range(T):
        D=delta[:,t,refs]-mu[t,refs]
        C=D.T@(D*sw[:,None]);C=(C+C.T)/2
        ridge=c['cov_ridge']*max(np.trace(C)/R,1e-12)+1e-12+c['ref_noise']**2
        S=C+np.eye(R)*ridge;inv=np.linalg.pinv(S,rcond=1e-15);invR[t]=inv
        A=delta[:,t,:]-mu[t]
        cross=A.T@(D*sw[:,None]);G=cross@inv;gain[t]=G
        variance=np.sum(A*A*sw[:,None],axis=0)-np.sum(G*cross,axis=1)
        cv[t]=np.maximum(variance,0)+c['target_noise']**2
    return mu,gain,cv,invR


def predict(nom,loads,y,delta_basis,state_w,pbar,c):
    # delta_basis state x time x output x bogie, N candidate nom x time x output
    delta=np.einsum('stob,b->sto',delta_basis,pbar,optimize=True)
    mu,gain,cv,invR=weighted_discrepancy(delta,state_w,REF,c)
    N,T,O=nom.shape;ll=np.zeros(N);ti=np.arange(0,T,c['likelihood_stride'])
    for t in ti:
        inn=y[t,REF][None,:]-nom[:,t,REF]-mu[t,REF][None,:]
        ll -= .5*c['likelihood_power']*np.sum((inn@invR[t])*inn,axis=1)
    ll-=ll.max();w=np.exp(ll);w/=w.sum();ess=1/np.sum(w*w)
    mean=np.zeros((T,O));sd=np.zeros((T,O))
    for t in range(T):
        inn=y[t,REF][None,:]-nom[:,t,REF]-mu[t,REF][None,:]
        a=nom[:,t,:]+mu[t][None,:]+inn@gain[t].T
        mean[t]=w@a
        second=w@(a*a)
        sd[t]=np.sqrt(np.maximum(second-mean[t]**2+cv[t],0))
    return dict(weights=w,ess=float(ess),mean=mean,sd=sd,lo=mean-c['interval_z']*sd,hi=mean+c['interval_z']*sd,
                pbar=w@loads,mu_delta=mu,gain=gain,condvar=cv,delta=delta,likelihood_indices=ti)


def scores(y,mean,lo,hi):
    rows=[]
    for idx,name in zip(TARGET,TARGET_NAMES):
        obs=y[:,idx];pred=mean[:,idx];rng=np.ptp(obs)
        rmse=float(np.sqrt(np.mean((pred-obs)**2)));nrmse=rmse/rng if rng>0 else np.nan
        r=float(np.corrcoef(pred,obs)[0,1]) if np.std(pred)>0 and np.std(obs)>0 else np.nan
        cov=float(np.mean((obs>=lo[:,idx])&(obs<=hi[:,idx])))
        amp=float(np.ptp(pred)/rng) if rng>0 else np.nan
        bias=float(np.mean(pred-obs))
        rows.append(dict(target=name,rmse_microstrain=rmse,nrmse=nrmse,correlation=r,coverage=cov,predicted_to_measured_range_ratio=amp,mean_bias_microstrain=bias))
    return pd.DataFrame(rows)


def main():
    start=time.time();assert self_test()<1e-12
    fe=loadmat(ROOT/'results'/'python_model.mat',simplify_cells=True);cfg=json.loads(str(fe['config_json']))
    nodes=np.asarray(fe['nodes'],float);shells=np.asarray(fe['shells'],float);P=np.asarray(fe['modes'],float);freq0=np.asarray(fe['frequencies_Hz'],float).ravel()
    summary=json.loads((SROOT/'results'/'shaker_ensemble_summary.json').read_text());states=pd.read_csv(SROOT/'results'/'shaker_weighted_states_all.csv');base=summary['baseline']
    ib=summary['stabil_mode_mapping']['first_global_bending_fe_index']-1;it=summary['stabil_mode_mapping']['first_global_torsion_fe_index']-1
    sel=np.where(freq0<=25.0)[0];sel=np.unique(np.r_[sel,ib,it]);Psel=P[:,sel]
    Bstrain,audit=build_bridgewest_strain_modal(nodes,shells,P,cfg,sel)
    pd.DataFrame(audit).to_csv(OUT/'strain_observation_map_audit.csv',index=False)
    # measured passage
    z=loadmat(PROCESSED,simplify_cells=True);tobs=np.asarray(z['time'],float).ravel();y=np.asarray(z['Y'],float);speed=float(z['speed'])
    if list(np.asarray(z['channels']).astype(str)) != ['str_x_BA_uf_031','str_x_BB_uf_031','str_x_BA_uf_053','str_x_BB_uf_053','str_x_BA_uf_058','str_x_BA_uf_070','str_x_BB_uf_081']:
        raise ValueError('Unexpected processed channel order')
    # internal time aligns exactly with 50 Hz tobs
    dt=.0025;ti=np.arange(tobs[0],tobs[-1]+dt/2,dt);sample_idx=np.rint((tobs-tobs[0])/dt).astype(int)
    if np.max(np.abs(ti[sample_idx]-tobs))>1e-9:raise ValueError('Time grids do not align')
    ax=np.array([0,2.6,19,21.6,25.9,28.5,44.9,47.5,51.8,54.4,70.8,73.4,77.7,80.3,96.7,99.3],float)
    pm=force_history(nodes,shells,Psel,sel,ti,speed,np.asarray(cfg['bearing_x_m'],float),ax,float(np.mean(cfg['girder_y_m'])))
    # nominal + all shaker states
    om,zeta,ob,zb=modal_parameters(freq0,sel,ib,it,states,base)
    nominal_basis=simulate_basis(pm,Bstrain,ob[None,:],zb[None,:],dt,sample_idx)[0].astype(float)*1e6 # T x output x bogie, microstrain per N
    state_basis=simulate_basis(pm,Bstrain,om,zeta,dt,sample_idx).astype(float)*1e6
    delta_basis=state_basis-nominal_basis[None,:,:,:]
    # load-only candidate library: deliberately ignore old stiffness/factor columns
    dr=loadmat(DRAWS,simplify_cells=True);loads=np.asarray(dr['loads'],float)*1e3
    nom=np.einsum('tob,nb->nto',nominal_basis,loads,optimize=True)
    c=dict(cov_ridge=.05,ref_noise=.5,target_noise=.5,likelihood_stride=10,likelihood_power=1.0,interval_z=1.959964)
    # reference-only fixed-point update of the load vector used to scale the structural discrepancy covariance
    pbar=loads.mean(0);trace=[]
    for k in range(6):
        res=predict(nom,loads,y,delta_basis,states.weight.to_numpy(),pbar,c)
        new=res['pbar'];change=float(np.max(np.abs(new-pbar))/1e3)
        trace.append(dict(iteration=k+1,max_bogie_change_kN=change,ess=res['ess'],**{f'p{j+1}_kN':new[j]/1e3 for j in range(8)}))
        pbar=new
        if change<.05:break
    pd.DataFrame(trace).to_csv(OUT/'reference_only_iteration_trace.csv',index=False)
    sc=scores(y,res['mean'],res['lo'],res['hi']);sc.to_csv(OUT/'run40_target_scores.csv',index=False)
    # references diagnostics; target values are used only below, after posterior is frozen.
    refrows=[]
    for idx in REF:
        obs=y[:,idx];pred=res['mean'][:,idx];rng=np.ptp(obs)
        refrows.append(dict(channel=CHANNELS[idx],nrmse=float(np.sqrt(np.mean((pred-obs)**2))/rng),correlation=float(np.corrcoef(pred,obs)[0,1]),coverage=float(np.mean((obs>=res['lo'][:,idx])&(obs<=res['hi'][:,idx])))))
    pd.DataFrame(refrows).to_csv(OUT/'run40_reference_fit.csv',index=False)
    # save histories
    hist={'time_s':tobs}
    for i,name in enumerate(CHANNELS):
        hist[f'{name}_measured']=y[:,i];hist[f'{name}_mean']=res['mean'][:,i];hist[f'{name}_lo']=res['lo'][:,i];hist[f'{name}_hi']=res['hi'][:,i]
    pd.DataFrame(hist).to_csv(OUT/'run40_histories.csv',index=False)
    loadtab=pd.DataFrame({'bogie':np.arange(1,9),'prior_mean_kN':loads.mean(0)/1e3,'posterior_mean_kN':res['pbar']/1e3,'posterior_sd_kN':np.sqrt(res['weights']@((loads/1e3-res['pbar'][None,:]/1e3)**2))})
    loadtab.to_csv(OUT/'run40_posterior_loads.csv',index=False)
    # plot four protected targets
    fig,axs=plt.subplots(2,2,figsize=(11,7),constrained_layout=True);axs=axs.ravel()
    for axp,idx,name in zip(axs,TARGET,TARGET_NAMES):
        axp.fill_between(tobs,res['lo'][:,idx],res['hi'][:,idx],alpha=.2,label='Nominal 95% interval')
        axp.plot(tobs,y[:,idx],label='Measured withheld')
        axp.plot(tobs,res['mean'][:,idx],label='Shaker-informed prediction')
        row=sc.loc[sc.target==name].iloc[0]
        axp.set_title(f'{name}: NRMSE={row.nrmse:.3f}, r={row.correlation:.3f}')
        axp.set_xlabel('Time from bridge entry (s)');axp.set_ylabel('Strain (µε)');axp.grid(True,alpha=.2)
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,bbox_to_anchor=(.5,-.02))
    fig.savefig(OUT/'run40_shaker_informed_targets.png',dpi=220,bbox_inches='tight');plt.close(fig)
    # posterior load plot
    fig,axp=plt.subplots(figsize=(8.5,4.5),constrained_layout=True)
    x=np.arange(1,9);axp.errorbar(x,loadtab.posterior_mean_kN,yerr=loadtab.posterior_sd_kN,fmt='o',capsize=3,label='Posterior')
    axp.plot(x,loadtab.prior_mean_kN,'--',label='Candidate-library mean')
    axp.set_xlabel('Bogie');axp.set_ylabel('Load (kN)');axp.set_xticks(x);axp.grid(True,alpha=.2);axp.legend()
    fig.savefig(OUT/'run40_posterior_loads.png',dpi=220);plt.close(fig)
    out=dict(event='2025_run40',speed_m_s=speed,observations=len(tobs),nominal_load_candidates=len(loads),shaker_states=len(states),shaker_effective_state_count=float(1/np.sum(states.weight.to_numpy()**2)),selected_dynamic_modes=len(sel),posterior_ess=res['ess'],posterior_mean_bogie_loads_kN=(res['pbar']/1e3).tolist(),iteration_count=len(trace),final_max_bogie_change_kN=trace[-1]['max_bogie_change_kN'],targets=sc.to_dict(orient='records'),references=refrows,elapsed_s=time.time()-start,
             caveat='Run 40 is a development/diagnostic passage previously inspected; these scores are not a sealed prospective validation. Target histories were opened only after the reference-conditioned posterior was frozen in this execution.')
    (OUT/'run40_summary.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
