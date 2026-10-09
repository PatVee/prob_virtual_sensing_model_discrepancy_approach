"""Numerical translation of the supplied Stabil shell4/DKT formulation.

Source routines and author comments are retained in dependencies/Stabil.
Rectangular, planar shells only. Beam elements are Euler-Bernoulli with
consistent translational mass and torsional inertia, no rotary bending inertia.
"""
import numpy as np
from scipy.linalg import block_diag
from q_dkt import q_dkt

def shape4(xi,eta):
    N=.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)])
    d=.25*np.array([[-(1-eta),1-eta,1+eta,-(1+eta)],[-(1-xi),-(1+xi),1+xi,1-xi]])
    return N,d

def dkt(nodes,h,E,nu,rho):
    b=nodes[[1,2,0],1]-nodes[[2,0,1],1]
    c=nodes[[2,0,1],0]-nodes[[1,2,0],0]
    det=b[1]*c[2]-b[2]*c[1]
    assert det>0
    D=E/(12*(1-nu**2))*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
    XX=h**3*np.array([[420,140,140],[140,70,35],[140,35,70]])
    Q=q_dkt(b,c,det)
    K=Q.T@np.kron(D,XX)/(det*840)@Q
    M=np.zeros((9,9))
    for L in np.array([[4,1,1],[1,4,1],[1,1,4]])/6:
        a=L[[1,2,0]];z=L[[2,0,1]];prod=np.prod(L)
        N=np.zeros(9)
        N[0::3]=L+L*L*a+L*L*z-L*a*a-L*z*z
        N[1::3]=b[[1,2,0]]*(z*L*L+prod/2)-b[[2,0,1]]*(L*L*a+prod/2)
        N[2::3]=c[[1,2,0]]*(z*L*L+prod/2)-c[[2,0,1]]*(L*L*a+prod/2)
        M+=np.outer(N,N)*rho*h*det/6
    return K,M

def shell_frame(nodes):
    ex=nodes[1]-nodes[0];ex/=np.linalg.norm(ex)
    ez=np.cross(nodes[2]-nodes[0],nodes[3]-nodes[1]);ez/=np.linalg.norm(ez)
    ey=np.cross(ez,ex);ex=np.cross(ey,ez)
    r=np.array([ex/np.linalg.norm(ex),ey/np.linalg.norm(ey),ez])
    loc=(nodes-nodes.mean(axis=0))@r.T
    assert np.max(np.abs(loc[:,2]))<1e-9
    return r,loc

def shell4(nodes,h,E,nu,rho):
    r,loc=shell_frame(nodes)
    T=block_diag(*([r]*8));K=np.zeros((24,24));M=K.copy()
    ids=np.array([0,1,6,7,12,13,18,19])
    D=E*h/(1-nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
    J0=.25*np.array([[-1,1,1,-1],[-1,-1,1,1]])@loc[:,:2]
    k1=np.zeros((8,8));k12=np.zeros((8,4));k2=np.zeros((4,4))
    for xi,eta in np.array([[-1,-1],[1,-1],[1,1],[-1,1]])/np.sqrt(3):
        N,d=shape4(xi,eta);J=d@loc[:,:2];det=np.linalg.det(J)
        dn=np.linalg.solve(J,d);B=np.zeros((3,8))
        B[0,0::2]=dn[0];B[1,1::2]=dn[1];B[2,0::2]=dn[1];B[2,1::2]=dn[0]
        de=np.linalg.solve(J0,np.linalg.det(J0)/det*np.diag([-2*xi,-2*eta]))
        BX=np.zeros((3,4));BX[0,[0,2]]=de[0];BX[1,[1,3]]=de[1];BX[2,[0,2]]=de[1];BX[2,[1,3]]=de[0]
        k1+=B.T@D@B*det;k12+=B.T@D@BX*det;k2+=BX.T@D@BX*det
        NN=np.zeros((2,8));NN[0,0::2]=N;NN[1,1::2]=N
        M[np.ix_(ids,ids)]+=NN.T@NN*rho*h*det
    K[np.ix_(ids,ids)]=k1-k12@np.linalg.solve(k2,k12.T)
    for seq in [[0,1,2],[0,2,3],[0,1,3],[1,2,3]]:
        kt,mt=dkt(loc[seq],h,E,nu,rho)
        ii=(6*np.array(seq)[:,None]+np.array([2,3,4])).ravel()
        K[np.ix_(ii,ii)]+=kt/2;M[np.ix_(ii,ii)]+=mt/2
    diag=np.diag(K);kd=np.min(diag[diag>0])/1000
    for j in [5,11,17,23]:K[j,j]=kd
    return T.T@K@T,T.T@M@T

def beam(nodes,sec,E,nu,rho):
    A,J,Iy,Iz=sec;delta=nodes[1]-nodes[0];L=np.linalg.norm(delta);ex=delta/L
    # Reference direction gives local z global-up for horizontal beams.
    ref=np.array([0.,0.,1.])
    if abs(ex@ref)>.9:ref=np.array([0.,1.,0.])
    ey=np.cross(ref,ex);ey/=np.linalg.norm(ey);ez=np.cross(ex,ey)
    T=block_diag(*([np.array([ex,ey,ez])]*4))
    K=np.zeros((12,12));M=K.copy()
    for ids,k,m in [([0,6],E*A/L,rho*A*L),([3,9],E/(2*(1+nu))*J/L,rho*(Iy+Iz)*L)]:
        K[np.ix_(ids,ids)]+=k*np.array([[1,-1],[-1,1]])
        M[np.ix_(ids,ids)]+=m/6*np.array([[2,1],[1,2]])
    k=np.array([[12,6*L,-12,6*L],[6*L,4*L*L,-6*L,2*L*L],[-12,-6*L,12,-6*L],[6*L,2*L*L,-6*L,4*L*L]])
    m=np.array([[156,22*L,54,-13*L],[22*L,4*L*L,13*L,-3*L*L],[54,13*L,156,-22*L],[-13*L,-3*L*L,-22*L,4*L*L]])
    for ids,I,sg in [([1,5,7,11],Iz,np.ones(4)),([2,4,8,10],Iy,np.array([1,-1,1,-1]))]:
        K[np.ix_(ids,ids)]+=E*I/L**3*k*np.outer(sg,sg)
        M[np.ix_(ids,ids)]+=rho*A*L/420*m*np.outer(sg,sg)
    return T.T@K@T,T.T@M@T
