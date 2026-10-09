"""Python counterpart of the Stabil se_shell4 observation path used by Bridge West.

This is not a replacement constitutive model. It is a direct translation of the
Stabil routines bundled with the package, used so the run-40 proof can be
executed in the delivery environment where MATLAB is unavailable.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import block_diag

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'python'))
from q_dkt import q_dkt


def _shape4(xi, eta):
    N = .25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)])
    dxi = .25*np.array([-(1-eta),1-eta,1+eta,-(1+eta)])
    deta = .25*np.array([-(1-xi),-(1+xi),1+xi,1-xi])
    return N, dxi, deta


def _frame(nodes):
    nx = nodes[1]-nodes[0]; nx /= np.linalg.norm(nx)
    v24 = nodes[3]-nodes[1]; v13 = nodes[2]-nodes[0]
    nz = np.cross(v13,v24); nz /= np.linalg.norm(nz)
    ny = np.cross(nz,nx); nx = np.cross(ny,nz)
    t = np.array([nx/np.linalg.norm(nx),ny/np.linalg.norm(ny),nz])
    local = (nodes-nodes.mean(0))@t.T
    if np.max(np.abs(local[:,2])) > 1e-8:
        raise ValueError('Bridge West observation translation assumes coplanar shell elements.')
    return t, local


def _vtrans_stress(t):
    th=np.zeros((6,6))
    for i in range(3):
        th[i,:3]=t[i,:]**2
        th[i,3]=2*t[i,0]*t[i,1]
        th[i,4]=2*t[i,1]*t[i,2]
        th[i,5]=2*t[i,2]*t[i,0]
    for row,(i,j) in enumerate([(0,1),(1,2),(2,0)],3):
        th[row,:3]=t[i,:]*t[j,:]
        th[row,3]=t[i,0]*t[j,1]+t[i,1]*t[j,0]
        th[row,4]=t[i,1]*t[j,2]+t[i,2]*t[j,1]
        th[row,5]=t[i,2]*t[j,0]+t[i,0]*t[j,2]
    return th


def _selcs_shell4(local,h,E,nu,U):
    gp=np.array([[-1,-1],[1,-1],[1,1],[-1,1]],float)/np.sqrt(3)
    D=E/(1-nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
    Dm=D*h
    J0=np.vstack([.25*np.array([-1,1,1,-1])@local[:,:2],
                   .25*np.array([-1,-1,1,1])@local[:,:2]])
    K21=np.zeros((4,8));K2=np.zeros((4,4))
    for xi,eta in gp:
        _,dxi,deta=_shape4(xi,eta)
        J=np.array([[dxi@local[:,0],dxi@local[:,1]],
                    [deta@local[:,0],deta@local[:,1]]])
        dNi=np.linalg.solve(J,np.vstack([dxi,deta]))
        Bm=np.zeros((3,8));Bm[0,0::2]=dNi[0];Bm[1,1::2]=dNi[1]
        Bm[2,0::2]=dNi[1];Bm[2,1::2]=dNi[0]
        Nex_der=np.array([[-2*xi,0],[0,-2*eta]])
        dNex=np.linalg.solve(J0,(np.linalg.det(J0)/np.linalg.det(J))*Nex_der.T)
        Bex=np.zeros((3,4));Bex[0,[0,2]]=dNex[0];Bex[1,[1,3]]=dNex[1]
        Bex[2,[0,2]]=dNex[1];Bex[2,[1,3]]=dNex[0]
        K21 += Bex.T@Dm@Bm*np.linalg.det(J)
        K2 += Bex.T@Dm@Bex*np.linalg.det(J)
    mem=np.array([0,1,6,7,12,13,18,19])
    Uex=-np.linalg.solve(K2,K21)@U[mem]
    Sg=np.zeros((12,U.shape[1]))
    for ig,(xi,eta) in enumerate(gp):
        _,dxi,deta=_shape4(xi,eta)
        J=np.array([[dxi@local[:,0],dxi@local[:,1]],
                    [deta@local[:,0],deta@local[:,1]]])
        dNi=np.linalg.solve(J,np.vstack([dxi,deta]))
        Bm=np.zeros((3,8));Bm[0,0::2]=dNi[0];Bm[1,1::2]=dNi[1]
        Bm[2,0::2]=dNi[1];Bm[2,1::2]=dNi[0]
        Nex_der=np.array([[-2*xi,0],[0,-2*eta]])
        dNex=np.linalg.solve(J0,(np.linalg.det(J0)/np.linalg.det(J))*Nex_der.T)
        Bex=np.zeros((3,4));Bex[0,[0,2]]=dNex[0];Bex[1,[1,3]]=dNex[1]
        Bex[2,[0,2]]=dNex[1];Bex[2,[1,3]]=dNex[0]
        Sg[3*ig:3*ig+3]=D@Bm@U[mem]+D@Bex@Uex
    g=1/abs(gp[0,0])/2;I=np.eye(3)
    extrap=np.block([
        [(1+g)*I,-.5*I,(1-g)*I,-.5*I],
        [-.5*I,(1+g)*I,-.5*I,(1-g)*I],
        [(1-g)*I,-.5*I,(1+g)*I,-.5*I],
        [-.5*I,(1-g)*I,-.5*I,(1+g)*I]])
    Sm=extrap@Sg
    Sb=np.zeros((40,U.shape[1]));seq=np.array([[0,1,2],[2,3,0],[3,0,1],[1,2,3]])
    L=np.array([[1,0,0],[1,1,0],[1,0,1]],float)
    for sq in seq:
        nod=local[sq]
        ind=np.r_[6*sq[0]+np.arange(2,5),6*sq[1]+np.arange(2,5),6*sq[2]+np.arange(2,5)]
        Ut=U[ind]
        b=np.array([nod[1,1]-nod[2,1],nod[2,1]-nod[0,1],nod[0,1]-nod[1,1]])
        c=np.array([nod[2,0]-nod[1,0],nod[0,0]-nod[2,0],nod[1,0]-nod[0,0]])
        det=b[1]*c[2]-b[2]*c[1];Q=q_dkt(b,c,det)
        chi1=L@Q[0:3]@Ut/det;chi2=L@Q[3:6]@Ut/det;chi12=L@Q[6:9]@Ut/det
        nd=sq[1]
        for iz,zeta in enumerate([1,-1]):
            vx=h*h/12*(D[0,0]*b[1:3]@Q[1:3]+D[0,1]*b[1:3]@Q[4:6]+D[2,2]*c[1:3]@Q[7:9])@Ut/det**2
            vy=h*h/12*(D[1,0]*c[1:3]@Q[1:3]+D[1,1]*c[1:3]@Q[4:6]+D[2,2]*b[1:3]@Q[7:9])@Ut/det**2
            start=5*nd+20*iz
            Sb[start:start+3]=(zeta*h/2)*D@np.vstack([chi1[1],chi2[1],chi12[1]])
            Sb[start+3]=vy;Sb[start+4]=vx
    S=np.zeros((72,U.shape[1]))
    for node in range(4):
        S[6*node+0]=Sm[3*node+0]+Sb[5*node+0]
        S[24+6*node+0]=Sm[3*node+0]
        S[48+6*node+0]=Sm[3*node+0]+Sb[20+5*node+0]
        S[6*node+1]=Sm[3*node+1]+Sb[5*node+1]
        S[24+6*node+1]=Sm[3*node+1]
        S[48+6*node+1]=Sm[3*node+1]+Sb[20+5*node+1]
        S[6*node+3]=Sm[3*node+2]+Sb[5*node+2]
        S[24+6*node+3]=Sm[3*node+2]
        S[48+6*node+3]=Sm[3*node+2]+Sb[20+5*node+2]
        S[24+6*node+4]=1.5*Sb[5*node+3]
        S[24+6*node+5]=1.5*Sb[5*node+4]
    return S


def shell_stress_operator(nodes,h,E,nu):
    t,local=_frame(nodes)
    T=block_diag(*([t]*8))
    Sl=_selcs_shell4(local,h,E,nu,T)  # U_global = I; U_local = T
    theta=_vtrans_stress(t.T)
    Sg=np.zeros_like(Sl)
    for block in range(12):
        Sg[6*block:6*block+6]=theta@Sl[6*block:6*block+6]
    return Sg


def build_bridgewest_strain_modal(nodes,shells,P,config,mode_indices):
    """Replicate bw_lorenzen_observation_maps.m and return 7 x nMode."""
    names=['A031','B031','A053_lower','B053','A058','A070','B081']
    xd=np.array([5.62,5.62,9.72,9.72,10.74,12.8,14.85])
    girder=np.array([0,1,0,1,0,0,1])
    E=config['E_Pa'];nu=config['nu'];L=config['length_m'];ys=np.asarray(config['girder_y_m']);bf=config['flange_width_m']
    Psel=P[:,mode_indices]
    Bm=np.zeros((7,len(mode_indices)))
    audit=[]
    for k in range(7):
        x=L-xd[k];g=girder[k];signSide=2*(g+1)-3
        y=ys[g]+signSide*bf/4
        group=3;offset=48
        if k==2:group=2;offset=0
        candidates=[]
        for e in np.flatnonzero(shells[:,5].astype(int)==group):
            ids=shells[e,:4].astype(int);p=nodes[ids]
            if x>=p[:,0].min()-1e-8 and x<=p[:,0].max()+1e-8 and y>=p[:,1].min()-1e-8 and y<=p[:,1].max()+1e-8:
                candidates.append(e)
        if not candidates:raise ValueError(f'{names[k]} outside flange mesh')
        row=np.zeros(len(mode_indices))
        for e in candidates:
            ids=shells[e,:4].astype(int);p=nodes[ids];h=shells[e,4]
            xi=2*(x-p[:,0].min())/(p[:,0].max()-p[:,0].min())-1
            eta=2*(y-p[:,1].min())/(p[:,1].max()-p[:,1].min())-1
            N=.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)])
            Sop=shell_stress_operator(p,h,E,nu)
            estrain=np.zeros((4,24))
            for j in range(4):
                ii=offset+6*j
                estrain[j]=(Sop[ii]-nu*(Sop[ii+1]+Sop[ii+2]))/E
            local_row=N@estrain
            dofs=(6*ids[:,None]+np.arange(6)).ravel()
            row += local_row@Psel[dofs]/len(candidates)
        Bm[k]=row
        audit.append(dict(channel=names[k],x_model_m=float(x),y_model_m=float(y),shell_group=group,surface_offset=offset,candidate_count=len(candidates)))
    return Bm, audit


def self_test():
    E=210e9;nu=.3;h=.02;eps=1e-5
    nodes=np.array([[0,0,0],[1,0,0],[1,.5,0],[0,.5,0]],float)
    U=np.zeros(24)
    for i,p in enumerate(nodes):U[6*i]=eps*p[0];U[6*i+1]=-nu*eps*p[1]
    S=shell_stress_operator(nodes,h,E,nu)
    stress=S@U
    vals=[]
    for off in [0,24,48]:
        for j in range(4):
            ii=off+6*j;vals.append((stress[ii]-nu*(stress[ii+1]+stress[ii+2]))/E)
    err=max(abs(np.asarray(vals)-eps))
    if err>1e-12:raise AssertionError(err)
    return err

if __name__=='__main__':
    print('self_test_error',self_test())
