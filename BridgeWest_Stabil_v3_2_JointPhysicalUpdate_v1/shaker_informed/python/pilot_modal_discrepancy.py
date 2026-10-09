from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[2]
SROOT=ROOT/'shaker_informed'
fe=loadmat(ROOT/'results'/'python_model.mat',simplify_cells=True)
c=json.loads(str(fe['config_json']))
nodes=np.asarray(fe['nodes'],float); shells=np.asarray(fe['shells'],float); P=np.asarray(fe['modes'],float); freq0=np.asarray(fe['frequencies_Hz'],float).ravel()
top=np.asarray(fe['top'],int); xx=np.asarray(fe['x'],float).ravel()
# scipy loadmat preserves the Python-exported zero-based node IDs.
sc=pd.read_json(SROOT/'results'/'shaker_ensemble_summary.json',typ='series')
states=pd.read_csv(SROOT/'results'/'shaker_weighted_states_all.csv')
summary=json.loads((SROOT/'results'/'shaker_ensemble_summary.json').read_text())
base=summary['baseline']; ib=summary['stabil_mode_mapping']['first_global_bending_fe_index']-1; it=summary['stabil_mode_mapping']['first_global_torsion_fe_index']-1
sel=np.where(freq0<=25.0)[0]; sel=np.unique(np.r_[sel,ib,it]); nm=len(sel)
# Observation: vertical displacement at the bearing-span midpoint on both upper-flange centrelines.
xmid=np.mean(c['bearing_x_m'])
def interp_mode(g):
    j=np.searchsorted(xx,xmid)-1;j=max(0,min(j,len(xx)-2));a=(xmid-xx[j])/(xx[j+1]-xx[j])
    ids=top[g,[j,j+1]]
    return (1-a)*P[6*ids[0]+2,sel]+a*P[6*ids[1]+2,sel]
obs=np.vstack([interp_mode(0),interp_mode(1)]) # 2 x nm

# Moving-load modal force history; 8 unit-bogie bases, track resultant at centreline.
v=160/3.6; dt=0.005; pre=0.5; post=2.0; a,b=c['bearing_x_m']; span=b-a
ax=np.array([0,2.6,19,21.6,25.9,28.5,44.9,47.5,51.8,54.4,70.8,73.4,77.7,80.3,96.7,99.3],float)
t=np.arange(-pre,(span+ax.max())/v+post+dt/2,dt); ytrack=np.mean(c['girder_y_m'])
deck=np.where(shells[:,5].astype(int)==4)[0]
# Cache modal vertical ordinates at all nodes.
phi_z=P[2::6,:][:,sel]
pm=np.zeros((len(t),nm,8))
def add_point(target,xp,yp,force):
    pick=None;best=1e99
    for ee in deck:
        ids=shells[ee,:4].astype(int); p=nodes[ids]; xmin,xmax=p[:,0].min(),p[:,0].max();ymin,ymax=p[:,1].min(),p[:,1].max()
        if xp>=xmin-1e-10 and xp<=xmax+1e-10 and yp>=ymin-1e-10 and yp<=ymax+1e-10:
            d=(xp-p[:,0].mean())**2+(yp-p[:,1].mean())**2
            if d<best:best=d;pick=ee
    if pick is None: return
    ids=shells[pick,:4].astype(int);p=nodes[ids];xmin,xmax=p[:,0].min(),p[:,0].max();ymin,ymax=p[:,1].min(),p[:,1].max()
    xi=2*(xp-xmin)/(xmax-xmin)-1;eta=2*(yp-ymin)/(ymax-ymin)-1
    N=.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)])
    target += force*(N[:,None]*phi_z[ids]).sum(0)
for i,tt in enumerate(t):
    xlead=a+v*tt
    for ia,off in enumerate(ax):
        x=xlead-off
        if a<=x<=b:add_point(pm[i,:,ia//2],x,ytrack,-0.5)

def modal_state(f1,f2,z1,z2):
    freq=freq0.copy();freq[ib]=f1;freq[it]=f2;om=2*np.pi*freq
    A=np.array([[1/(2*om[ib]),om[ib]/2],[1/(2*om[it]),om[it]/2]])
    alpha,beta=np.linalg.solve(A,np.array([z1,z2]));zeta=alpha/(2*om)+beta*om/2;zeta=np.clip(zeta,.002,.08);zeta[ib]=z1;zeta[it]=z2
    return om[sel],zeta[sel]
def simulate(f1,f2,z1,z2):
    om,z=modal_state(f1,f2,z1,z2);kk=om**2;cc=2*z*om
    beta=.25;gamma=.5;a0=1/(beta*dt**2);a1=gamma/(beta*dt);a2=1/(beta*dt);a3=1/(2*beta)-1;a4=gamma/beta-1;a5=dt*(gamma/(2*beta)-1);den=kk+a0+a1*cc
    q=np.zeros((nm,8));qd=q.copy();qdd=q.copy(); y=np.zeros((len(t),2,8))
    for j in range(len(t)):
        force=pm[j]
        if j==0:qdd=force-cc[:,None]*qd-kk[:,None]*q
        else:
            pe=force+a0*q+a2*qd+a3*qdd+cc[:,None]*(a1*q+a4*qd+a5*qdd)
            qn=pe/den[:,None];an=a0*(qn-q)-a2*qd-a3*qdd;vn=qd+dt*((1-gamma)*qdd+gamma*an);q,qd,qdd=qn,vn,an
        y[j]=obs@q
    return y
nom=simulate(base['mode1_frequency_Hz'],base['mode2_frequency_Hz'],base['mode1_damping_ratio'],base['mode2_damping_ratio'])
load=np.full(8,270e3); ynom=np.tensordot(nom,load,axes=(2,0))*1e3 # mm
# Evaluate all 84 states.
peak=[];rms=[];deltas=[]
for r in states.itertuples(index=False):
    y=simulate(r.f1_Hz,r.f2_Hz,r.zeta1,r.zeta2); yy=np.tensordot(y,load,axes=(2,0))*1e3; d=yy-ynom;deltas.append(d)
    peak.append(np.max(np.abs(d),axis=0));rms.append(np.sqrt(np.mean(d*d,axis=0)))
deltas=np.array(deltas);peak=np.array(peak);rms=np.array(rms);w=states.weight.to_numpy();
# weighted response discrepancy SD through time, then summarize.
mu=np.tensordot(w,deltas,axes=(0,0));var=np.tensordot(w,(deltas-mu)**2,axes=(0,0));sd=np.sqrt(var)
out={
 'pilot_only':'Python mechanics check of modal discrepancy using midspan BA/BB displacement; MATLAB Stabil shell-strain recovery remains the production path.',
 'speed_kmh':160.0,'mean_bogie_load_kN':270.0,'time_step_s':dt,'selected_modes':int(nm),
 'nominal_baseline':base,
 'weighted_mean_peak_abs_discrepancy_mm_BA_BB':(w[:,None]*peak).sum(0).tolist(),
 'weighted_mean_rms_discrepancy_mm_BA_BB':(w[:,None]*rms).sum(0).tolist(),
 'maximum_weighted_pointwise_sd_mm_BA_BB':sd.max(axis=0).tolist(),
 'nominal_peak_abs_displacement_mm_BA_BB':np.max(np.abs(ynom),axis=0).tolist(),
 'finite':bool(np.isfinite(deltas).all()),'nonzero':bool(np.max(np.abs(deltas))>0),
}
(SROOT/'results'/'python_dynamic_pilot_summary.json').write_text(json.dumps(out,indent=2))
# save compact time history of weighted discrepancy mean/sd
pd.DataFrame({'time_s':t,'mean_BA_mm':mu[:,0],'sd_BA_mm':sd[:,0],'mean_BB_mm':mu[:,1],'sd_BB_mm':sd[:,1]}).to_csv(SROOT/'results'/'python_dynamic_pilot_discrepancy.csv',index=False)
print(json.dumps(out,indent=2))
