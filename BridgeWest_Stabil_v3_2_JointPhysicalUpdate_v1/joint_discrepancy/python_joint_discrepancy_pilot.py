from __future__ import annotations
from pathlib import Path
import json,sys,time,math
import numpy as np, pandas as pd
from scipy.io import loadmat
from scipy.optimize import minimize
from scipy.special import ndtr,ndtri,log_ndtr,ndtri_exp

ROOT=Path('/mnt/data/bridgewest_compare/v3pkg/BridgeWest_Stabil_v3_2_MassAudit_Plots/BridgeWest_Stabil_v3_2_MassAudit_Plots')
V3=Path('/mnt/data/bridgewest_compare/v3res/BridgeWest_Stabil_Shaker_v3_results')
OUT=Path('/mnt/data/bridgewest_joint_discrepancy_pilot'); OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'shaker_informed/python'))
from stabil_shell_observation import build_bridgewest_strain_modal

DCFG=json.loads((ROOT/'shaker_informed/discrepancy_config.json').read_text())
PCFG=json.loads((ROOT/'production/config/prediction_config.json').read_text())
STATES=pd.read_csv(ROOT/'shaker_informed/results/shaker_weighted_states_all.csv')
SUMMARY=json.loads((ROOT/'shaker_informed/results/shaker_ensemble_summary.json').read_text())
BASE=SUMMARY['baseline']; REF=np.array([0,3],int); TARGET=np.array([4,5,1,6],int)
TARGET_NAMES=['A058','A070','B031','B081']; AXLES=np.asarray(DCFG['axle_offsets_m'],float)
GAUGE=1.435
# Epistemic load-transfer states. Equal weights; none is selected using target/load/displacement validation.
TRANSFER_STATES=[
    ('two_rail_point',0.0,0.0),
    ('spread_short',0.6,0.30),
    ('spread_medium',1.2,0.45),
    ('spread_long',2.4,0.60),
]
TRANSFER_W=np.ones(len(TRANSFER_STATES))/len(TRANSFER_STATES)
# Zero-mean (in arithmetic response-gain space) residual amplitude layers.
COMPLIANCE_GAINS=np.array([0.80,1.00,1.20]); COMPLIANCE_W=np.array([0.25,0.50,0.25])
STRAIN_RECOVERY_GAINS=np.array([0.92,1.00,1.08]); STRAIN_RECOVERY_W=np.array([0.25,0.50,0.25])


def displacement_modal(fe,sel):
    top=np.asarray(fe['top']).astype(int); xnodes=np.asarray(fe['x'],float).ravel(); P=np.asarray(fe['modes'],float)
    cfg=json.loads(str(fe['config_json'])); a,b=cfg['bearing_x_m']; xs=np.linspace(a,b,41)
    out=[]
    for g in (0,1):
        for xx in xs:
            k=np.searchsorted(xnodes,xx,side='right')-1; k=max(0,min(k,len(xnodes)-2))
            aa=(xx-xnodes[k])/(xnodes[k+1]-xnodes[k]); ids=top[g,[k,k+1]].astype(int)
            out.append((1-aa)*P[6*ids[0]+2,sel]+aa*P[6*ids[1]+2,sel])
    # dataset direction used by package plotter
    return np.asarray(out), np.asarray(cfg['length_m']-xs)


def modal_parameters(freq0,sel,ib,it,states,base):
    n=len(states); freq=np.tile(freq0[sel],(n,1)); pb=int(np.flatnonzero(sel==ib)[0]); pt=int(np.flatnonzero(sel==it)[0])
    freq[:,pb]=states.f1_Hz.to_numpy(); freq[:,pt]=states.f2_Hz.to_numpy(); om=2*np.pi*freq; z=np.empty_like(om)
    for s,r in enumerate(states.itertuples(index=False)):
        o1=2*np.pi*r.f1_Hz; o2=2*np.pi*r.f2_Hz
        A=np.array([[1/(2*o1),o1/2],[1/(2*o2),o2/2]])
        al,be=np.linalg.solve(A,np.array([r.zeta1,r.zeta2])); z[s]=np.clip(al/(2*om[s])+be*om[s]/2,.002,.08); z[s,pb]=r.zeta1; z[s,pt]=r.zeta2
    fb=freq0[sel].copy(); fb[pb]=base['mode1_frequency_Hz']; fb[pt]=base['mode2_frequency_Hz']; ob=2*np.pi*fb
    o1=2*np.pi*base['mode1_frequency_Hz']; o2=2*np.pi*base['mode2_frequency_Hz']; A=np.array([[1/(2*o1),o1/2],[1/(2*o2),o2/2]])
    al,be=np.linalg.solve(A,np.array([base['mode1_damping_ratio'],base['mode2_damping_ratio']])); zb=np.clip(al/(2*ob)+be*ob/2,.002,.08); zb[pb]=base['mode1_damping_ratio']; zb[pt]=base['mode2_damping_ratio']
    return om,z,ob,zb


def deck_grid(nodes,shells,Psel,bearing):
    de=np.flatnonzero(shells[:,5].astype(int)==4); ids=np.unique(shells[de,:4].astype(int)); x=nodes[ids,0]; y=nodes[ids,1]
    ux=np.unique(np.round(x,10)); uy=np.unique(np.round(y,10)); xb={};yb={}
    for i,v in enumerate(ux): xb[round(float(v),10)]=((ux[i-1]+v)/2 if i else v,(v+ux[i+1])/2 if i<len(ux)-1 else v)
    for i,v in enumerate(uy): yb[round(float(v),10)]=((uy[i-1]+v)/2 if i else v,(v+uy[i+1])/2 if i<len(uy)-1 else v)
    xL=[];xR=[];yB=[];yT=[]
    for xx,yy in zip(x,y):
        a,c=xb[round(float(xx),10)]; b,d=yb[round(float(yy),10)]; xL.append(a);xR.append(c);yB.append(b);yT.append(d)
    return dict(ids=ids,xL=np.array(xL),xR=np.array(xR),yB=np.array(yB),yT=np.array(yT),phi=Psel[6*ids+2],bearing=np.array(bearing,float),deck=de)


def point_force(grid,nodes,shells,Psel,xp,yp):
    for e in grid['deck']:
        ids=shells[e,:4].astype(int); p=nodes[ids]
        if p[:,0].min()-1e-10<=xp<=p[:,0].max()+1e-10 and p[:,1].min()-1e-10<=yp<=p[:,1].max()+1e-10:
            xmin,xmax=p[:,0].min(),p[:,0].max(); ymin,ymax=p[:,1].min(),p[:,1].max(); xi=2*(xp-xmin)/(xmax-xmin)-1; eta=2*(yp-ymin)/(ymax-ymin)-1
            N=.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)])
            return -(N[:,None]*Psel[6*ids+2]).sum(0)
    return np.zeros(Psel.shape[1])


def area_force(grid,xp,yp,hx,hy):
    a,b=grid['bearing']; xl=max(xp-hx,a); xr=min(xp+hx,b)
    if xr<=xl:return np.zeros(grid['phi'].shape[1])
    dx=np.maximum(0,np.minimum(grid['xR'],xr)-np.maximum(grid['xL'],xl)); dy=np.maximum(0,np.minimum(grid['yT'],yp+hy)-np.maximum(grid['yB'],yp-hy)); area=dx*dy; s=area.sum()
    return -(area[:,None]/s*grid['phi']).sum(0) if s>0 else np.zeros(grid['phi'].shape[1])


def force_history(state,grid,nodes,shells,Psel,times,speed,bearing,ytrack):
    name,hx,hy=state; pm=np.zeros((len(times),Psel.shape[1],8),dtype=float); a,b=bearing
    for it,t in enumerate(times):
        xlead=a+speed*t
        for ia,off in enumerate(AXLES):
            xp=xlead-off
            if xp<a-1e-12 or xp>b+1e-12: continue
            bog=ia//2
            for yr in (ytrack-GAUGE/2,ytrack+GAUGE/2):
                if hx<=0: f=point_force(grid,nodes,shells,Psel,xp,yr)
                else: f=area_force(grid,xp,yr,hx,hy)
                pm[it,:,bog] += .25*f
    return pm


def simulate(pm,obs,om,zeta,dt,sidx):
    ns,nm=om.shape; nb=8; no=obs.shape[0]; kk=om**2; cc=2*zeta*om; be=.25;ga=.5
    a0=1/(be*dt**2);a1=ga/(be*dt);a2=1/(be*dt);a3=1/(2*be)-1;a4=ga/be-1;a5=dt*(ga/(2*be)-1); den=kk+a0+a1*cc
    q=np.zeros((ns,nm,nb));qd=np.zeros_like(q);qdd=np.zeros_like(q);out=np.empty((ns,len(sidx),no,nb),np.float32); mp={int(k):j for j,k in enumerate(sidx)}
    for j in range(len(pm)):
        f=pm[j][None]
        if j==0:qdd=f-cc[:,:,None]*qd-kk[:,:,None]*q
        else:
            pe=f+a0*q+a2*qd+a3*qdd+cc[:,:,None]*(a1*q+a4*qd+a5*qdd);qn=pe/den[:,:,None];an=a0*(qn-q)-a2*qd-a3*qdd;vn=qd+dt*((1-ga)*qdd+ga*an);q,qd,qdd=qn,vn,an
        if j in mp: out[:,mp[j]]=np.einsum('om,smb->sob',obs,q,optimize=True)
    return out


def component_stats(delta,weights):
    # delta: state x time x output
    w=np.asarray(weights,float);w=w/w.sum(); mu=np.einsum('s,sto->to',w,delta,optimize=True)
    return mu,delta-mu[None]


def joint_moments(component_sets,refs,ref_noise,target_noise_by_group,ridge_frac=.05):
    # component_sets: list of dicts {'groups':[delta_strain,delta_ba,delta_bb], 'weights':w}
    T=component_sets[0]['groups'][0].shape[1]; R=len(refs); means=[None]*3; centered=[]
    for k in range(3):
        means[k]=sum(component_stats(c['groups'][k],c['weights'])[0] for c in component_sets)
    # store centered arrays + weights per component/group once
    compc=[]
    for c in component_sets:
        gcent=[]
        for k in range(3):
            mu,cen=component_stats(c['groups'][k],c['weights']);gcent.append(cen)
        compc.append((gcent,np.asarray(c['weights'],float)/np.sum(c['weights'])))
    gains=[np.zeros((T,component_sets[0]['groups'][k].shape[2],R)) for k in range(3)]
    cvars=[np.zeros((T,component_sets[0]['groups'][k].shape[2])) for k in range(3)]
    invR=np.zeros((T,R,R)); Crr_all=np.zeros((T,R,R))
    for t in range(T):
        Crr=np.zeros((R,R))
        for gcent,w in compc:
            D=gcent[0][:,t,refs]; Crr += D.T@(D*w[:,None])
        Crr=(Crr+Crr.T)/2; Crr_all[t]=Crr
        ridge=ridge_frac*max(np.trace(Crr)/R,1e-12)+1e-12+ref_noise**2
        S=Crr+np.eye(R)*ridge; inv=np.linalg.pinv(S,rcond=1e-15); invR[t]=inv
        for k in range(3):
            O=gains[k].shape[1]; cross=np.zeros((O,R)); var=np.zeros(O)
            for gcent,w in compc:
                A=gcent[k][:,t,:]; D=gcent[0][:,t,refs]
                cross += A.T@(D*w[:,None]); var += np.sum(A*A*w[:,None],axis=0)
            G=cross@inv; gains[k][t]=G; cvars[k][t]=np.maximum(var-np.sum(G*cross,axis=1),0)+target_noise_by_group[k]**2
    return dict(mean=means,gain=gains,cvar=cvars,invR=invR,Crr=Crr_all)


def rtr(a,b,rng):
    if a>0:return -rtr(-b,-a,rng)
    if b>0:
        Fa=ndtr(a);Fb=ndtr(b);return ndtri(np.clip(Fa+(Fb-Fa)*rng.random(),1e-300,1-1e-16))
    la=log_ndtr(a);lb=log_ndtr(b);rr=np.exp(la-lb) if np.isfinite(la) else 0.;return ndtri_exp(lb+np.log(rr+(1-rr)*rng.random()))

def bmap(H,h,lo,hi):
    r=minimize(lambda x:.5*x@H@x-h@x,(lo+hi)/2,jac=lambda x:H@x-h,bounds=list(zip(lo,hi)),method='L-BFGS-B',options={'ftol':1e-12,'gtol':1e-8,'maxiter':1000});return r.x

def tg(H,h,lo,hi,seed,burn,keep):
    H=(H+H.T)/2; ev=np.linalg.eigvalsh(H); rg=max(0,1e-12*max(ev.max(),1)-ev.min()); H=H+np.eye(len(lo))*rg
    m=np.linalg.solve(H,h);x=bmap(H,h,lo,hi);rng=np.random.default_rng(seed);out=np.empty((keep,len(lo)));n=0
    for s in range(burn+keep):
        for i in range(len(lo)):
            sd=np.sqrt(1/H[i,i]); cm=m[i]-(H[i]@(x-m)-H[i,i]*(x[i]-m[i]))/H[i,i]; x[i]=cm+sd*rtr((lo[i]-cm)/sd,(hi[i]-cm)/sd,rng)
        if s>=burn:out[n]=x;n+=1
    return out.mean(0),np.cov(out,rowvar=False,ddof=1),out,m

def posterior(Bref,y,mu,inv,lo,hi,seed,burn,keep):
    H=np.zeros((8,8));h=np.zeros(8)
    for t in range(0,len(y),PCFG['likelihood_time_stride']):
        B=Bref[t];S=inv[t];d=y[t]-mu[t];H+=B.T@S@B;h+=B.T@S@d
    return tg(H,h,lo,hi,seed,burn,keep)


def make_components(pbarN,B0,BDYN,BTRANS):
    # B*: groups; B0[k] T x O x 8; BDYN[k] S x T x O x 8; BTRANS[k] L x T x O x 8
    comps=[]
    dyn=[np.einsum('stob,b->sto',BDYN[k]-B0[k][None],pbarN,optimize=True) for k in range(3)]
    comps.append(dict(name='shaker_dynamic',groups=dyn,weights=STATES.weight.to_numpy(float)))
    tr=[np.einsum('ltob,b->lto',BTRANS[k]-B0[k][None],pbarN,optimize=True) for k in range(3)]
    comps.append(dict(name='load_transfer',groups=tr,weights=TRANSFER_W))
    cg=[]
    for k in range(3):
        base=np.einsum('tob,b->to',B0[k],pbarN,optimize=True)
        cg.append(np.stack([(g-1)*base for g in COMPLIANCE_GAINS],axis=0))
    comps.append(dict(name='common_compliance',groups=cg,weights=COMPLIANCE_W))
    rg=[]
    base_s=np.einsum('tob,b->to',B0[0],pbarN,optimize=True)
    rg.append(np.stack([(g-1)*base_s for g in STRAIN_RECOVERY_GAINS],axis=0))
    for k in (1,2): rg.append(np.zeros((len(STRAIN_RECOVERY_GAINS),B0[k].shape[0],B0[k].shape[1])))
    comps.append(dict(name='strain_recovery',groups=rg,weights=STRAIN_RECOVERY_W))
    return comps


def predict_groups(B0,mom,q,C,yref):
    out=[]; Br=B0[0][:,REF,:]
    for k in range(3):
        T,O,_=B0[k].shape; mn=np.zeros((T,O));sd=np.zeros((T,O))
        for t in range(T):
            B=B0[k][t]*1000;Brr=Br[t]*1000;G=mom['gain'][k][t];A=B-G@Brr
            cc=mom['mean'][k][t]+G@(yref[t]-mom['mean'][0][t,REF])
            mn[t]=cc+A@q; v=np.einsum('oi,ij,oj->o',A,C,A,optimize=True);sd[t]=np.sqrt(np.maximum(v+mom['cvar'][k][t],0))
        out.append((mn,sd))
    return out


def load_inputs(event):
    z=loadmat(V3/'permitted_inputs'/f'{event}.mat',simplify_cells=True)
    return np.asarray(z['time_s'],float).ravel(),np.asarray(z['reference_strain_microstrain'],float),float(z['speed_m_s'])

def static_bogies(year):
    w=pd.read_csv(ROOT/'audit/verified_static_bogie_loads.csv');return w[w.campaign==year].sort_values('bogie').static_load_kN.to_numpy(float)


def event_bases(event):
    t,yref,v=load_inputs(event);dt=DCFG['simulation_dt_s'];ti=np.arange(t[0],t[-1]+dt/2,dt);sidx=np.rint((t-t[0])/dt).astype(int);ytrack=float(np.mean(CFG['girder_y_m']))
    pms=[]
    for state in TRANSFER_STATES:pms.append(force_history(state,GRID,NODES,SHELLS,PSEL,ti,v,np.asarray(CFG['bearing_x_m'],float),ytrack))
    pmbar=sum(w*pm for w,pm in zip(TRANSFER_W,pms))
    B0=[];BDYN=[];BTRANS=[]
    for obs,scale in [(BSTRAIN,1e6),(BDISP[:41],1e3),(BDISP[41:],1e3)]:
        b0=simulate(pmbar,obs,OB[None],ZB[None],dt,sidx)[0].astype(float)*scale
        bd=simulate(pmbar,obs,OM,ZETA,dt,sidx).astype(float)*scale
        bt=np.stack([simulate(pm,obs,OB[None],ZB[None],dt,sidx)[0].astype(float)*scale for pm in pms],axis=0)
        B0.append(b0);BDYN.append(bd);BTRANS.append(bt)
    return t,yref,v,B0,BDYN,BTRANS


def infer(event,full=False):
    t,yref,v,B0,BDYN,BTRANS=event_bases(event)
    lo=np.full(8,PCFG['load_lower_kN']);hi=np.full(8,PCFG['load_upper_kN']);p=(lo+hi)/2
    bi=PCFG['gibbs_burn_iteration'] if full else 300;ki=PCFG['gibbs_keep_iteration'] if full else 1500
    bf=PCFG['gibbs_burn_final'] if full else 600;kf=PCFG['gibbs_keep_final'] if full else 4000
    trace=[]
    for itn in range(PCFG['fixed_point_max_iterations']):
        comps=make_components(p*1000,B0,BDYN,BTRANS);mom=joint_moments(comps,REF,PCFG['reference_noise_sigma_microstrain'],[PCFG['target_noise_sigma_microstrain'],PCFG['displacement_noise_sigma_mm'],PCFG['displacement_noise_sigma_mm']],PCFG['cov_ridge'])
        raw,C,samp,uncon=posterior(B0[0][:,REF,:]*1000,yref,mom['mean'][0][:,REF],mom['invR'],lo,hi,PCFG['gibbs_seed']+itn,bi,ki)
        pp=(1-PCFG['fixed_point_relaxation'])*p+PCFG['fixed_point_relaxation']*raw;ch=float(np.max(abs(pp-p)));trace.append((itn+1,ch,float(pp.sum())));p=pp
        if ch<PCFG['fixed_point_tolerance_kN']:break
    for ipass in range(2):
        comps=make_components(p*1000,B0,BDYN,BTRANS);mom=joint_moments(comps,REF,PCFG['reference_noise_sigma_microstrain'],[PCFG['target_noise_sigma_microstrain'],PCFG['displacement_noise_sigma_mm'],PCFG['displacement_noise_sigma_mm']],PCFG['cov_ridge'])
        q,C,samp,uncon=posterior(B0[0][:,REF,:]*1000,yref,mom['mean'][0][:,REF],mom['invR'],lo,hi,PCFG['gibbs_seed']+100+ipass,bf,kf);p=q
    # freeze before validation
    pred=predict_groups(B0,mom,q,C,yref)
    total_samples=samp.sum(axis=1);qtot=np.quantile(total_samples,[.025,.5,.975]);yr=int(event[:4]);qs=static_bogies(yr)
    result=dict(event=event,speed_m_s=v,q=q,C=C,samples=samp,total=float(q.sum()),total_q025=float(qtot[0]),total_q50=float(qtot[1]),total_q975=float(qtot[2]),static_total=float(qs.sum()),trace=trace,pred=pred,t=t,B0=B0,mom=mom)
    # validation only from here
    z=loadmat(V3/'validation_targets'/f'{event}.mat',simplify_cells=True);yt=np.asarray(z['Y'],float);flags=np.asarray(z['flags']).astype(str)
    srows=[];smn,ssd=pred[0]
    for j,name in enumerate(TARGET_NAMES):
        yy=yt[:,j];pp=smn[:,TARGET[j]];rng=np.ptp(yy);loi=pp-1.95996398454*ssd[:,TARGET[j]];hii=pp+1.95996398454*ssd[:,TARGET[j]]
        srows.append(dict(target=name,flag=flags[j],nrmse=float(np.sqrt(np.mean((pp-yy)**2))/rng),correlation=float(np.corrcoef(pp,yy)[0,1]),range_ratio=float(np.ptp(pp)/rng),coverage=float(np.mean((yy>=loi)&(yy<=hii)))))
    result['strain_scores']=srows
    snap=V3/'displacement_validation'/f'{event}_measured_PLS_snapshot.csv';drows=[]
    if snap.exists():
        truth=pd.read_csv(snap)
        for gi,g in enumerate(['BA','BB']):
            dmn,dsd=pred[1+gi];x=XDISP.copy();mid=(x.min()+x.max())/2;jm=np.argmin(abs(x-mid));it=np.argmax(abs(dmn[:,jm]));pr=dmn[it];sd=dsd[it];ord=np.argsort(x);xp=x[ord];pr=pr[ord];sd=sd[ord]
            qq=truth[truth.axis==g].sort_values('x_m');xx=qq.x_m.to_numpy(float);yy=qq.measured_mm.to_numpy(float);pi=np.interp(xx,xp,pr);si=np.interp(xx,xp,sd);rng=np.ptp(yy);loi=pi-1.95996398454*si;hii=pi+1.95996398454*si
            drows.append(dict(girder=g,peak_ratio=float(np.max(abs(pi))/np.max(abs(yy))),spatial_nrmse=float(np.sqrt(np.mean((pi-yy)**2))/rng),correlation=float(np.corrcoef(pi,yy)[0,1]),profile_coverage=float(np.mean((yy>=loi)&(yy<=hii))),pred_peak_mm=float(np.min(pi)),meas_peak_mm=float(np.min(yy)),pred_peak_lo_mm=float(np.min(loi)),pred_peak_hi_mm=float(np.min(hii))))
    result['disp_scores']=drows
    return result

# model setup
FE=loadmat(ROOT/'results/python_model.mat',simplify_cells=True);CFG=json.loads(str(FE['config_json']));NODES=np.asarray(FE['nodes'],float);SHELLS=np.asarray(FE['shells'],float);P=np.asarray(FE['modes'],float);FREQ=np.asarray(FE['frequencies_Hz'],float).ravel()
ib=SUMMARY['stabil_mode_mapping']['first_global_bending_fe_index']-1;it=SUMMARY['stabil_mode_mapping']['first_global_torsion_fe_index']-1;SEL=np.where(FREQ<=DCFG['dynamic_max_frequency_Hz'])[0];SEL=np.unique(np.r_[SEL,ib,it]);PSEL=P[:,SEL]
BSTRAIN,_=build_bridgewest_strain_modal(NODES,SHELLS,P,CFG,SEL);BDISP,XDISP=displacement_modal(FE,SEL);OM,ZETA,OB,ZB=modal_parameters(FREQ,SEL,ib,it,STATES,BASE);GRID=deck_grid(NODES,SHELLS,PSEL,CFG['bearing_x_m'])

if __name__=='__main__':
    events=sys.argv[1:] or ['2022_run20','2025_run22','2025_run40']
    rows=[];srows=[];drows=[]
    meta=dict(transfer_states=[dict(name=n,hx_m=hx,hy_m=hy,weight=float(w)) for (n,hx,hy),w in zip(TRANSFER_STATES,TRANSFER_W)],compliance_gains=COMPLIANCE_GAINS.tolist(),compliance_weights=COMPLIANCE_W.tolist(),strain_recovery_gains=STRAIN_RECOVERY_GAINS.tolist(),strain_recovery_weights=STRAIN_RECOVERY_W.tolist(),shaker_states=len(STATES),protocol='A031/B053 inference only; target strain/static weighing/PLS validation after posterior freeze')
    (OUT/'joint_discrepancy_config.json').write_text(json.dumps(meta,indent=2))
    for ev in events:
        st=time.time();r=infer(ev,full=False);elapsed=time.time()-st
        print('DONE',ev,'total',r['total'],'95',r['total_q025'],r['total_q975'],'static',r['static_total'],'sec',elapsed,flush=True)
        rows.append(dict(event=ev,total_kN=r['total'],total_q025_kN=r['total_q025'],total_q975_kN=r['total_q975'],static_total_kN=r['static_total'],total_error_pct=100*(r['total']/r['static_total']-1),static_inside_95=bool(r['total_q025']<=r['static_total']<=r['total_q975']),n_lower=int(np.sum(r['q']<=182)),n_upper=int(np.sum(r['q']>=358)),elapsed_s=elapsed,**{f'q{i+1}_kN':r['q'][i] for i in range(8)}))
        for z in r['strain_scores']:srows.append(dict(event=ev,**z))
        for z in r['disp_scores']:drows.append(dict(event=ev,**z))
    pd.DataFrame(rows).to_csv(OUT/'pilot_loads.csv',index=False);pd.DataFrame(srows).to_csv(OUT/'pilot_strain_scores.csv',index=False);pd.DataFrame(drows).to_csv(OUT/'pilot_displacement_scores.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(drows).to_string(index=False))
