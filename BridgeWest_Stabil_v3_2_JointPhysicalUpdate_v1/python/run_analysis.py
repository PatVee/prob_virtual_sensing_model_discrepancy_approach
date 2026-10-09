"""Rebuild, check mechanics and compare modal families without traffic fitting."""
import csv,json,copy,time
import numpy as np
from scipy.io import savemat
from scipy.sparse.linalg import spsolve
from model import *

def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def compare(spectrum):
    targets=list(csv.DictReader((ROOT/'evidence/thesis_modal_targets.csv').open()))
    families=['bending','torsion','bending','torsion'];lobes=[1,1,2,2];used=set()
    for i,r in enumerate(targets):
        r.update(fe_index='',fe_frequency_Hz='',frequency_error_percent='',assignment_status='unassigned')
        if i==4:
            r['assignment_status']='undefined in thesis; no automatic assignment';continue
        cand=[s for s in spectrum if s['family']==families[i] and s['lobes']==lobes[i] and s['fe_index'] not in used]
        if cand:
            t=float(r['frequency_Hz']);s=min(cand,key=lambda s:abs(np.log(s['frequency_Hz']/t)))
            used.add(s['fe_index']);r.update(fe_index=s['fe_index'],fe_frequency_Hz=s['frequency_Hz'],frequency_error_percent=100*(s['frequency_Hz']/t-1),assignment_status='family-based candidate; measured MAC unavailable')
    return targets

def checks(m):
    K=m['K'];M=m['M'];n=K.shape[0];P=m['P'];nmode=P.shape[1]
    symmetryK=np.sqrt((K-K.T).power(2).sum()/K.power(2).sum())
    symmetryM=np.sqrt((M-M.T).power(2).sum()/M.power(2).sum())
    total=[]
    for j in range(3):
        q=np.zeros(n);q[j::6]=1;total.append(float(q@M@q))
    target=m['c']['length_m']*m['c']['mass_per_length_total_kg_m']
    # Equal downward deck force to check reactions and bounded deflection.
    p=np.zeros(n);d=np.flatnonzero(m['added_mass']>0)
    p[6*d+2]=-1e3/len(d)
    q=np.zeros(n);f=m['free'];q[f]=spsolve(K[f][:,f],p[f]);reaction=K@q-p
    force_balance=float(np.sum(reaction[2::6])-np.dot(m['foundation'][2::6],q[2::6])+np.sum(p[2::6]))
    orth=float(np.linalg.norm(P.T@(M@P)-np.eye(nmode)))
    r=dict(nodes=len(m['nodes']),shell_elements=len(m['shells']),beam_elements=len(m['beams']),free_dofs=len(f),
        stiffness_symmetry_error=float(symmetryK),mass_symmetry_error=float(symmetryM),
        translation_mass_kg=total,nominal_total_mass_kg=target,steel_mass_kg=m['steel_mass'],
        steel_density_kg_m3=m['rho'],residual_track_ballast_mass_kg=m['residual_mass'],
        walkway_surface_mass_kg=m['walkway_mass'],mass_orthogonality_error=orth,
        max_eigenpair_residual=float(max(m['residual'])),minimum_frequency_Hz=float(m['freq'][0]),
        highest_computed_frequency_Hz=float(m['freq'][-1]),bearing_displacement_under_1kN_m=float(max(abs(q[m['fixed']]))),
        downward_static_peak_m=float(np.min(q[2::6])),vertical_force_balance_N=force_balance,
        foundation_total_stiffness_N_m=float(np.sum(m['foundation'])),geometry_verified=False,
        MATLAB_executed=False,execution='Python translation of supplied Stabil element formulation')
    assert symmetryK<1e-12 and symmetryM<1e-12 and max(abs(np.array(total)-target))/target<1e-9
    assert orth<1e-6 and max(m['residual'])<1e-5 and abs(force_balance)<1e-5
    return r

def shaker_compare(m):
    rows=list(csv.DictReader((ROOT/'evidence/shaker_shapes_2025.csv').open()));out=[];spec=modal_classification(m)
    for k,kind in enumerate(['bending','torsion'],1):
        candidate=[r for r in spec if r['family']==kind and r['lobes']==1]
        if not candidate:continue
        index=min(candidate,key=lambda r:r['frequency_Hz'])['fe_index']-1
        a=[];b=[]
        for row in rows:
            g=0 if row['girder']=='A' else 1
            x=1+float(row['x_FE_m'])
            a.append(np.interp(x,m['x'],m['P'][6*m['top'][g]+2,index]));b.append(float(row[f'mode{k}_STW']))
        a=np.array(a);b=np.array(b);mac=float((a@b)**2/((a@a)*(b@b)))
        out.append(dict(mode=kind,fe_index=index+1,MAC_2025_shaker=mac,role='different campaign from thesis; diagnostic only'))
    return out

def main():
    c=json.loads((ROOT/'config/model.json').read_text());start=time.time()
    m=solve(assemble(mesh(c)));s=modal_classification(m);r=checks(m)
    write_csv(ROOT/'results/modal_spectrum.csv',s);write_csv(ROOT/'results/thesis_comparison.csv',compare(s))
    write_csv(ROOT/'results/shaker_shape_check.csv',shaker_compare(m))
    export_mesh(m,ROOT/'results/python_model.mat')
    (ROOT/'results/mechanics_checks.json').write_text(json.dumps(r,indent=2))
    print('Nominal',compare(s),flush=True)
    variants=[]
    for name,settings in [('no_upper_flange_restraint',{'top_foundation_N_m3':0}),
                          ('longitudinal_upper_flange_restraint',{'top_spring_direction':'x'}),
                          ('half_length_patch',{'top_patch_length_m':.5}),
                          ('refined_mesh',{'max_dx_m':c['max_dx_m']/2,'max_dy_m':c['max_dy_m']/2,'max_dz_m':c['max_dz_m']/2})]:
        cc={**c,**settings};mm=solve(assemble(mesh(cc)));sp=modal_classification(mm)
        for row in compare(sp):variants.append({'scenario':name,**row})
        print(name,[(row['type'],row['fe_frequency_Hz']) for row in compare(sp)],flush=True)
        if name=='refined_mesh':
            # Matched low global modes at common top-centreline coordinates.
            export_mesh(mm,ROOT/'results/refined_model.mat')
    write_csv(ROOT/'results/assumption_and_mesh_sensitivity.csv',variants)
    print('Completed in',round(time.time()-start,1),'s',flush=True)

if __name__=='__main__':main()
