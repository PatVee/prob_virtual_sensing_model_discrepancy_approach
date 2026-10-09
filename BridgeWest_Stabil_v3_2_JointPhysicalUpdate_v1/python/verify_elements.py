"""Independent analytical checks of shell/beam mechanics used in the pilot."""
import json
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from model import ROOT
from elements import shell4,beam

E=2.1e11;nu=.3;rho=7850.;h=.02
p=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]])
K,M=shell4(p,h,E,nu,rho)
q=np.zeros((4,6));q[:,2]=.5*p[:,0]**2;q[:,4]=-p[:,0]
energy=q.ravel()@K@q.ravel();expected=E*h**3/(12*(1-nu**2))
assert abs(energy/expected-1)<1e-10

def plate(n):
    a=2.;b=1.4;nodes=np.array([[x,y,0] for x in np.linspace(0,a,n+1) for y in np.linspace(0,b,n+1)])
    rows=[];cols=[];ks=[];ms=[]
    for i in range(n):
        for j in range(n):
            ids=np.array([i*(n+1)+j,(i+1)*(n+1)+j,(i+1)*(n+1)+j+1,i*(n+1)+j+1])
            ke,me=shell4(nodes[ids],h,E,nu,rho);ii=(6*ids[:,None]+np.arange(6)).ravel()
            rows.extend(np.repeat(ii,24));cols.extend(np.tile(ii,24));ks.extend(ke.ravel());ms.extend(me.ravel())
    nd=6*len(nodes);k=coo_matrix((ks,(rows,cols)),shape=(nd,nd)).tocsr();m=coo_matrix((ms,(rows,cols)),shape=(nd,nd)).tocsr()
    fixed=list((6*np.arange(len(nodes))[:,None]+np.array([0,1,5])).ravel())
    edge=np.flatnonzero((nodes[:,0]==0)|(nodes[:,0]==a)|(nodes[:,1]==0)|(nodes[:,1]==b));fixed.extend(6*edge+2)
    free=np.setdiff1d(np.arange(nd),fixed)
    vals=eigsh(k[free][:,free],k=3,M=m[free][:,free],sigma=0,return_eigenvectors=False)
    f=np.sqrt(min(vals))/(2*np.pi);analytic=np.pi/2*np.sqrt(expected/(rho*h))*(1/a**2+1/b**2)
    return dict(mesh=n,FE_Hz=float(f),analytical_Hz=float(analytic),error_percent=float(100*(f/analytic-1)))

plates=[plate(4),plate(8),plate(12)]
assert abs(plates[-1]['error_percent'])<2.5
assert abs(plates[-1]['error_percent'])<abs(plates[0]['error_percent'])
# Slender cantilever bending frequency: a different element and boundary case.
A=.01;I=1e-5;L=5.;n=12;nd=6*(n+1);K=np.zeros((nd,nd));M=K.copy()
for i in range(n):
    p=np.array([[L*i/n,0.,0.],[L*(i+1)/n,0.,0.]])
    k,m=beam(p,[A,2*I,I,I],E,nu,rho);ids=np.arange(6*i,6*i+12);K[np.ix_(ids,ids)]+=k;M[np.ix_(ids,ids)]+=m
v=eigh(K[6:,6:],M[6:,6:],eigvals_only=True);f=np.sqrt(v[0])/(2*np.pi)
analytic=1.8751040687**2/(2*np.pi*L**2)*np.sqrt(E*I/(rho*A))
assert abs(f/analytic-1)<1e-4
out=dict(shell_constant_curvature_energy_relative_error=float(energy/expected-1),simply_supported_plate=plates,
         cantilever_beam_FE_Hz=float(f),cantilever_analytical_Hz=float(analytic),status='passed',
         limitation='Python element translation checked analytically; not MATLAB runtime parity.')
(ROOT/'results/element_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
