"""Parameterised, provisional trough-bridge mesh and sparse FE assembly."""
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import eigsh
from scipy.io import savemat
from elements import shell4,beam,shape4,shell_frame

ROOT=Path(__file__).resolve().parents[1]

def subdivide(vals,step):
    vals=np.unique(np.round(vals,9));out=[]
    for a,b in zip(vals[:-1],vals[1:]):out.extend(np.linspace(a,b,int(np.ceil((b-a)/step-1e-8))+1)[:-1])
    return np.r_[out,vals[-1]]

def mesh(c):
    nodes=[];lookup={};shells=[];beams=[]
    def node(p):
        key=tuple(np.round(p,8))
        if key not in lookup:lookup[key]=len(nodes);nodes.append(key)
        return lookup[key]
    L=c['length_m'];a,b=c['bearing_x_m'];ys=c['girder_y_m'];h=c['girder_depth_m'];bf=c['flange_width_m'];zd=c['deck_z_m']
    cross=np.unique(np.r_[0,np.linspace(a,b,int(round((b-a)/c['crossbeam_spacing_m']))+1),L])
    patch=c['top_patch_length_m']
    xx=subdivide([0,patch,L-patch,L,a,b,*cross],c['max_dx_m'])
    zz=subdivide([0,zd,h],c['max_dz_m']);ydeck=subdivide(ys,c['max_dy_m'])
    def surface(us,vs,point,t,group):
        grid=np.array([[node(point(u,v)) for v in vs] for u in us])
        for i in range(len(us)-1):
            for j in range(len(vs)-1):shells.append([*grid[[i,i+1,i+1,i],[j,j,j+1,j+1]],t,group])
        return grid
    for y in ys:
        surface(xx,zz,lambda x,z:[x,y,z],c['web_thickness_m'],1)
        for z,group in [(0,2),(h,3)]:
            surface(xx,[y-bf/2,y,y+bf/2],lambda x,Y:[x,Y,z],c['flange_thickness_m'],group)
    surface(xx,ydeck,lambda x,y:[x,y,zd],c['deck_thickness_m'],4)
    def section_i():
        d=c['crossbeam_depth_m'];w=c['crossbeam_width_m'];tw=c['crossbeam_web_m'];tf=c['crossbeam_flange_m']
        A=2*w*tf+(d-2*tf)*tw
        Iy=(w*d**3-(w-tw)*(d-2*tf)**3)/12
        Iz=2*tf*w**3/12+(d-2*tf)*tw**3/12
        J=(2*w*tf**3+(d-2*tf)*tw**3)/3
        return [A,J,Iy,Iz]
    sec=section_i()
    # Crossbeams tied at the deck surface. No unverified eccentric offsets.
    for x in cross:
        for y0,y1 in zip(ydeck[:-1],ydeck[1:]):beams.append([node([x,y0,zd]),node([x,y1,zd]),*sec,1])
    if c['include_walkway']:
        yout=ys[0]-c['walkway_width_m'];zw=zd
        for x in cross:
            beams.append([node([x,ys[0],zw]),node([x,yout,zw]),*sec,2])
        D=c['handrail_tube_diameter_m'];d=D-2*c['handrail_tube_wall_m']
        A=np.pi*(D**2-d**2)/4;I=np.pi*(D**4-d**4)/64;rail=[A,2*I,I,I]
        for x in cross:beams.append([node([x,yout,zw]),node([x,yout,zw+c['handrail_height_m']]),*rail,3])
        for x0,x1 in zip(cross[:-1],cross[1:]):
            beams.append([node([x0,yout,zw]),node([x1,yout,zw]),*sec,2])
            beams.append([node([x0,yout,zw+c['handrail_height_m']]),node([x1,yout,zw+c['handrail_height_m']]),*rail,3])
    nodes=np.array(nodes);shells=np.array(shells);beams=np.array(beams)
    fixed=[]
    for x in [a,b]:
        for y in ys:fixed.append(6*node([x,y,0])+2)
    fixed.extend([6*node([a,ys[0],0]),6*node([a,ys[0],0])+1,6*node([b,ys[0],0])+1])
    # Sensors lie on upper flange centreline; no arbitrary recovery depth.
    top=np.array([[node([x,y,h]) for x in xx] for y in ys])
    assert len(nodes)==len(lookup)
    return dict(nodes=nodes,shells=shells,beams=beams,fixed=np.unique(fixed),top=top,x=xx,c=c)

def assemble(m):
    c=m['c'];nodes=m['nodes'];n=6*len(nodes);E=c['E_Pa'];nu=c['nu'];rho=c['specific_weight_N_m3']/c['g_m_s2']
    rows=[];cols=[];ks=[];ms=[];cache={};steel_mass=0.;deck_areas=np.zeros(len(nodes));walk_areas=np.zeros(len(nodes));foundation=np.zeros(n)
    def add(ii,k,mm):
        rows.extend(np.repeat(ii,len(ii)));cols.extend(np.tile(ii,len(ii)));ks.extend(k.ravel());ms.extend(mm.ravel())
    for row in m['shells']:
        ids=row[:4].astype(int);p=nodes[ids];t=row[4];group=int(row[5]);ii=(6*ids[:,None]+np.arange(6)).ravel()
        key=(tuple(np.round((p-p[0]).ravel(),8)),t)
        if key not in cache:cache[key]=shell4(p,t,E,nu,rho)
        ke,me=cache[key];add(ii,ke,me)
        area=np.linalg.norm(np.cross(p[1]-p[0],p[3]-p[0]));steel_mass+=rho*t*area
        if group==4:deck_areas[ids]+=area/4
        if group==3 and (p[:,0].max()<=c['top_patch_length_m']+1e-7 or p[:,0].min()>=c['length_m']-c['top_patch_length_m']-1e-7):
            # Tributary area lumping of the distributed foundation; sum is k*A.
            axis={'x':0,'y':1,'z':2}[c['top_spring_direction']]
            foundation[6*ids+axis]+=c['top_foundation_N_m3']*area/4
    for row in m['beams']:
        ids=row[:2].astype(int);sec=row[2:6];p=nodes[ids];ii=(6*ids[:,None]+np.arange(6)).ravel()
        key=('beam',tuple(np.round((p-p[0]).ravel(),8)),tuple(sec))
        if key not in cache:cache[key]=beam(p,sec,E,nu,rho)
        ke,me=cache[key];add(ii,ke,me);steel_mass+=rho*sec[0]*np.linalg.norm(p[1]-p[0])
        if int(row[6])==2 and abs(p[1,0]-p[0,0])>1e-7:
            walk_areas[ids]+=np.linalg.norm(p[1]-p[0])*c['walkway_width_m']/2
    K=coo_matrix((ks,(rows,cols)),shape=(n,n)).tocsr()+diags(foundation)
    M=coo_matrix((ms,(rows,cols)),shape=(n,n)).tocsr()
    walkway_mass=walk_areas*c['walkway_surface_mass_kg_m2']
    residual_mass=c['mass_per_length_total_kg_m']*c['length_m']-steel_mass-walkway_mass.sum()
    if residual_mass<=0:raise ValueError('Steel and walkway exceed prescribed total mass')
    added=walkway_mass+residual_mass*deck_areas/deck_areas.sum()
    massdiag=np.zeros(n)
    for j in range(3):massdiag[j::6]=added
    M+=diags(massdiag)
    # Normal shell rotations can have zero physical mass. Keep K and use
    # shift-invert eigsh with semidefinite M; no artificial mass is introduced.
    free=np.setdiff1d(np.arange(n),m['fixed'])
    m.update(K=K,M=M,free=free,foundation=foundation,added_mass=added,steel_mass=steel_mass,
             residual_mass=residual_mass,walkway_mass=float(walkway_mass.sum()),rho=rho)
    return m

def solve(m,nm=None):
    f=m['free'];K=m['K'][f][:,f];M=m['M'][f][:,f];nm=nm or m['c']['number_of_modes']
    vals,V=eigsh(K,k=nm,M=M,sigma=-1e-4,which='LM',tol=2e-9,v0=np.sin(np.arange(len(f))+1))
    order=np.argsort(vals);vals=vals[order];V=V[:,order]
    if vals[0]<=0:raise ValueError(f'Nonpositive eigenvalue {vals[0]}')
    V/=np.sqrt(np.sum(V*(M@V),axis=0))
    full=np.zeros((6*len(m['nodes']),nm));full[f]=V
    residual=np.linalg.norm(K@V-(M@V)*vals,axis=0)/(np.linalg.norm(K@V,axis=0)+np.linalg.norm((M@V)*vals,axis=0))
    m.update(freq=np.sqrt(vals)/(2*np.pi),lam=vals,P=full,residual=residual)
    return m

def modal_classification(m):
    # Physical mode families from motion of the two upper-flange centrelines.
    x=m['x'];a,b=m['c']['bearing_x_m'];mask=(x>a+1e-6)&(x<b-1e-6)
    Z=m['P'][6*m['top']+2,:];records=[]
    for j,f in enumerate(m['freq']):
        A=Z[0,mask,j];B=Z[1,mask,j];energy=np.dot(A,A)+np.dot(B,B)
        sym=np.dot(A+B,A+B)/(2*energy+1e-30)
        anti=np.dot(A-B,A-B)/(2*energy+1e-30)
        transl=m['P'][:,j].reshape(-1,6)[:,:3]
        vertical=float(np.sum(transl[:,2]**2)/np.sum(transl**2))
        top_ratio=float(np.max(np.abs(np.r_[A,B]))/(np.max(np.abs(transl[:,2]))+1e-30))
        v=A+B if sym>=anti else A-B
        v=v[np.abs(v)>.08*np.max(np.abs(v))]
        lobes=1+np.sum(v[:-1]*v[1:]<0) if len(v) else 0
        family=('bending' if sym>=anti else 'torsion') if vertical>.5 and top_ratio>.25 else 'local/lateral'
        records.append(dict(fe_index=j+1,frequency_Hz=float(f),family=family,lobes=int(lobes),symmetric_fraction=float(sym),vertical_fraction=vertical,top_amplitude_fraction=top_ratio))
    return records

def export_mesh(m,path):
    savemat(path,dict(nodes=m['nodes'],shells=m['shells'],beams=m['beams'],fixed=m['fixed'],top=m['top'],x=m['x'],
                     foundation=m['foundation'],added_mass=m['added_mass'],frequencies_Hz=m['freq'],modes=m['P'],
                     K=m['K'],M=m['M'],config_json=json.dumps(m['c'])),do_compression=True)

if __name__=='__main__':
    c=json.loads((ROOT/'config/model.json').read_text());m=solve(assemble(mesh(c)))
    for row in modal_classification(m):print(row)
    export_mesh(m,ROOT/'results/python_model.mat')
