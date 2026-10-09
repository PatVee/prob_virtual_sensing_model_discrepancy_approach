from pathlib import Path
import sys,json,time
import numpy as np,pandas as pd
from scipy.io import loadmat
sys.path.insert(0,'/mnt/data/bridgewest_compare')
import joint_discrepancy_v3 as j
from stabil_shell_observation import build_bridgewest_strain_modal

OUT=Path('/mnt/data/bridgewest_joint_calibrated');OUT.mkdir(exist_ok=True)
# Production transfer mixture: distributed two-rail states only; point contact remains a diagnostic, not a production state.
j.TRANSFER_STATES=[('spread_short',0.6,0.30),('spread_medium',1.2,0.45),('spread_long',2.4,0.60)]
j.TRANSFER_W=np.ones(3)/3
# Independent 2024 single-shaker force-normalised modal-residue evidence.
OLD_RES1=2.142090694997057e-05;OLD_RES2=5.204858470581088e-05
OLD_TO_MEAS1=1.59272552937298;OLD_TO_MEAS2=0.945167979654921
MEAS1=OLD_RES1*OLD_TO_MEAS1;MEAS2=OLD_RES2*OLD_TO_MEAS2

def mode_at_top(P,g,xp,mode):
    top=np.asarray(j.FE['top']).astype(int);xn=np.asarray(j.FE['x'],float).ravel();k=np.searchsorted(xn,xp,side='right')-1;k=max(0,min(k,len(xn)-2));a=(xp-xn[k])/(xn[k+1]-xn[k]);ids=top[g,[k,k+1]];v=P[6*ids+2,mode];return float((1-a)*v[0]+a*v[1])
def residue(P,mode):
    left=j.CFG['bearing_x_m'][0];xs=[left+4.1,left+8.2];resp=np.array([mode_at_top(P,0,xs[0],mode),mode_at_top(P,1,xs[0],mode),mode_at_top(P,0,xs[1],mode),mode_at_top(P,1,xs[1],mode)]);force=mode_at_top(P,1,8.7,mode);return abs(force)*np.linalg.norm(resp)
R1=MEAS1/residue(j.P,j.ib);R2=MEAS2/residue(j.P,j.it)
P2=j.P.copy();P2[:,j.ib]*=np.sqrt(R1);P2[:,j.it]*=np.sqrt(R2)
# Patch joint module's force participation and observation maps consistently.
j.P=P2;j.PSEL=P2[:,j.SEL]
j.BSTRAIN,_=build_bridgewest_strain_modal(j.NODES,j.SHELLS,P2,j.CFG,j.SEL)
FE2=dict(j.FE);FE2['modes']=P2;j.BDISP,j.XDISP=j.displacement_modal(FE2,j.SEL);j.GRID=j.deck_grid(j.NODES,j.SHELLS,j.PSEL,j.CFG['bearing_x_m'])

# Fast 2022 strain-recovery calibration using static weighing and permitted A031/B053 only.
def nominal_ref_basis(event):
    t,y,v=j.load_inputs(event);dt=j.DCFG['simulation_dt_s'];ti=np.arange(t[0],t[-1]+dt/2,dt);sidx=np.rint((t-t[0])/dt).astype(int);ytrack=float(np.mean(j.CFG['girder_y_m']))
    pms=[j.force_history(st,j.GRID,j.NODES,j.SHELLS,j.PSEL,ti,v,np.asarray(j.CFG['bearing_x_m'],float),ytrack) for st in j.TRANSFER_STATES]
    pmbar=sum(w*pm for w,pm in zip(j.TRANSFER_W,pms));B=j.simulate(pmbar,j.BSTRAIN,j.OB[None],j.ZB[None],dt,sidx)[0].astype(float)*1e6
    return t,y,B

cal=[]
for run in [2,4,6,8,10,12,14,16,18,20,24,26]:
    ev=f'2022_run{run:02d}';t,y,B=nominal_ref_basis(ev);q=j.static_bogies(2022)*1000;pred=np.einsum('trb,b->tr',B[:,j.REF,:],q,optimize=True)
    ss=[]
    for k,name in enumerate(['A031','B053']):
        den=float(pred[:,k]@pred[:,k]);s=float(pred[:,k]@y[:,k]/den) if den>0 else np.nan;ss.append(s)
    pv=pred.reshape(-1);yv=y.reshape(-1);sc=float(pv@yv/(pv@pv))
    cal.append(dict(event=ev,scale_event=sc,scale_A031=ss[0],scale_B053=ss[1]))
caldf=pd.DataFrame(cal);SCALE=float(caldf.scale_event.median());caldf.to_csv(OUT/'2022_reference_scale_by_event.csv',index=False)
(OUT/'calibration.json').write_text(json.dumps(dict(modal_residue_scale_bending=R1,modal_residue_scale_torsion=R2,measured_residue_bending=MEAS1,measured_residue_torsion=MEAS2,strain_recovery_scale_2022_median=SCALE,scale_event_min=float(caldf.scale_event.min()),scale_event_max=float(caldf.scale_event.max()),policy='2024 shaker amplitude + 2022 static bogie weighing/A031-B053 only; 2025 static weighing, target strain and PLS not used'),indent=2))
print('R1,R2',R1,R2,'strain scale',SCALE,'range',caldf.scale_event.min(),caldf.scale_event.max(),flush=True)

orig_event_bases=j.event_bases
def scaled_event_bases(event):
    t,y,v,B0,BD,BT=orig_event_bases(event)
    B0[0]=B0[0]*SCALE;BD[0]=BD[0]*SCALE;BT[0]=BT[0]*SCALE
    return t,y,v,B0,BD,BT
j.event_bases=scaled_event_bases

# Execute 2025 validation after calibration is frozen.
events=sys.argv[1:] or ['2025_run22','2025_run40']
rows=[];sr=[];dr=[]
for ev in events:
    st=time.time();r=j.infer(ev,full=False);el=time.time()-st
    rows.append(dict(event=ev,total_kN=r['total'],total_q025_kN=r['total_q025'],total_q975_kN=r['total_q975'],static_total_kN=r['static_total'],total_error_pct=100*(r['total']/r['static_total']-1),static_inside_95=bool(r['total_q025']<=r['static_total']<=r['total_q975']),elapsed_s=el,**{f'q{i+1}_kN':r['q'][i] for i in range(8)}))
    for z in r['strain_scores']:sr.append(dict(event=ev,**z))
    for z in r['disp_scores']:dr.append(dict(event=ev,**z))
    print('DONE',ev,'total',r['total'],'95',r['total_q025'],r['total_q975'],'static',r['static_total'],'disp',r['disp_scores'],'sec',el,flush=True)
pd.DataFrame(rows).to_csv(OUT/'2025_load_validation.csv',index=False);pd.DataFrame(sr).to_csv(OUT/'2025_strain_validation.csv',index=False);pd.DataFrame(dr).to_csv(OUT/'2025_displacement_validation.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(dr).to_string(index=False))
