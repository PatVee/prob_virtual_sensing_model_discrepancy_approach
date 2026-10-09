"""Paper-style diagnostic figures. Reads results; does not run the solver."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

def save(fig,name):
    fig.savefig(ROOT/'figures'/f'{name}.png',dpi=190,bbox_inches='tight')
    fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)

def bridge(ax,s,mode=None):
    N=s['nodes'];F=s['shells'][:,:4].astype(int);B=s['beams'][:,:2].astype(int)
    if mode is None:
        P=N;colors=np.array(['#bac7cf','#647b8c','#647b8c','#b9bfc3'])[s['shells'][:,5].astype(int)-1]
    else:
        u=s['modes'][:,mode].reshape(-1,6)[:,:3];scale=.45/max(np.max(abs(u)),1e-20);P=N+scale*u
        z=u[F,2].mean(axis=1);colors=plt.cm.RdBu_r(.5+.5*z/(max(abs(z))+1e-30))
    az=np.deg2rad(-65);el=np.deg2rad(24)
    H=np.array([-np.sin(az),np.cos(az),0]);V=np.array([-np.sin(el)*np.cos(az),-np.sin(el)*np.sin(az),np.cos(el)])
    D=np.cross(H,V);proj=np.column_stack([P@H,P@V]);order=np.argsort(P[F].mean(axis=1)@D)
    coll=PolyCollection(proj[F[order]],facecolors=colors[order],edgecolors=(.15,.2,.25,.2),linewidths=.18,alpha=.86);ax.add_collection(coll)
    ax.add_collection(LineCollection(proj[B],colors='.25',linewidths=.5,alpha=.6))
    if mode is None:
        sup=np.unique(s['fixed'].ravel()//6);p=proj[sup];ax.scatter(*p.T,marker='^',s=42,color='#9a3e35',zorder=8)
        spring=np.unique(np.flatnonzero(s['foundation'].ravel())//6);p=proj[spring];ax.scatter(*p.T,s=11,color='#208375',zorder=9)
    ax.autoscale();ax.set_aspect('equal');ax.margins(.03);ax.axis('off')

def main():
    s=loadmat(ROOT/'results/python_model.mat');comp=list(csv.DictReader((ROOT/'results/thesis_comparison.csv').open()))
    fig=plt.figure(figsize=(14,5.3));ax=fig.add_axes([.03,.18,.94,.68]);bridge(ax,s)
    ax.set_title('Bridge West: provisional shell-and-beam reconstruction',pad=1,fontsize=16)
    fig.text(.08,.09,'Red triangles: bearing restraints     Green patches: upper-flange foundation springs',fontsize=11)
    fig.text(.08,.045,'Plate, stiffener and walkway dimensions remain assumptions; this is not the original SOFiSTiK mesh.',fontsize=10,color='.3')
    save(fig,'bw_lorenzen_geometry')
    fig=plt.figure(figsize=(15,8.4));fig.suptitle('Comparison with Lorenzen Table 7.2 — provisional geometry',fontsize=17,y=.98)
    for i,r in enumerate(comp):
        if r['fe_index']:
            ax=fig.add_subplot(3,2,i+1);j=int(r['fe_index'])-1;bridge(ax,s,j)
            qualifier='; mesh-sensitive candidate' if i==2 else '; candidate'
            ax.set_title(f"({chr(97+i)}) {r['type'].capitalize()}{qualifier}\nFE {float(r['fe_frequency_Hz']):.2f} Hz | measured {float(r['frequency_Hz']):.2f} Hz",fontsize=12,pad=0)
        else:
            ax=fig.add_subplot(3,2,i+1);ax.axis('off')
            ax.text(.07,.74,f"({chr(97+i)}) Thesis mode {i+1}: {r['type']}\nMeasured frequency {float(r['frequency_Hz']):.2f} Hz",fontsize=15,va='top')
            text=('No reliable family assignment\nin the provisional reconstruction.' if i==3 else 'Undefined in the thesis (MPC = 0.54).\nNo FE mode forced to match this frequency.')
            ax.text(.07,.37,text,fontsize=12,color='.3',va='top')
    ax=fig.add_subplot(3,2,6);ax.axis('off')
    ax.text(.07,.78,'Assessment',fontsize=15,weight='bold')
    ax.text(.07,.65,'First two shapes are stable under mesh refinement.\nHigher-mode mapping is unresolved.\n\nNo train passage or target response was fitted.\nThis model has not passed physical validation.',fontsize=12,va='top',linespacing=1.5)
    fig.subplots_adjust(top=.9,bottom=.06,hspace=.30,wspace=.14);save(fig,'bw_lorenzen_mode_comparison')
    fig,axs=plt.subplots(2,2,figsize=(12,7.5),constrained_layout=True)
    for i,(ax,r) in enumerate(zip(axs.flat,comp[:4])):
        if not r['fe_index']:
            ax.axis('off');ax.text(.1,.55,'Second torsion: no reliable assignment',transform=ax.transAxes,fontsize=12);continue
        j=int(r['fe_index'])-1;x=s['x'][0];A=s['modes'][6*s['top'][0]+2,j];B=s['modes'][6*s['top'][1]+2,j];scale=max(abs(np.r_[A,B]));A/=scale;B/=scale
        if A[np.argmax(abs(A))]<0:A=-A;B=-B
        ax.plot(x,A,'-',color='black',lw=1.7,label='BA');ax.plot(x,B,'--',color='#236a89',lw=1.7,label='BB')
        for xp in [1,17.4]:ax.axvline(xp,color='.6',lw=.8,ls=':')
        ax.axhline(0,color='.6',lw=.6);ax.set(xlim=(0,18.4),ylim=(-1.1,1.1),xlabel='Model coordinate x (m)',ylabel='Normalised vertical mode',title=f"({chr(97+i)}) {r['type'].capitalize()} — {float(r['fe_frequency_Hz']):.2f} Hz")
        ax.grid(alpha=.17);ax.legend(frameon=False,loc='best')
    save(fig,'bw_lorenzen_girder_profiles')
    spec=list(csv.DictReader((ROOT/'results/modal_spectrum.csv').open()))
    fig,ax=plt.subplots(figsize=(12,4),constrained_layout=True)
    for r in spec:
        f=float(r['frequency_Hz'])
        if f<=30:ax.plot([f,f],[.05,.42],color='#607d8b' if r['family']!='local/lateral' else '.78',lw=1.2)
    for i,r in enumerate(comp):
        f=float(r['frequency_Hz']);ax.plot([f,f],[.65,.9],color='black',lw=2)
        ax.text(f,.96,f'{i+1}\n{f:.2f}',ha='center',va='bottom',fontsize=11)
    ax.set(xlim=(0,30),ylim=(0,1.24),xlabel='Frequency (Hz)',yticks=[.24,.77],yticklabels=['Provisional FE spectrum','Thesis observations'],title='Observed modes cannot be assigned by the order of FE eigenvalues alone')
    ax.grid(axis='x',alpha=.18);save(fig,'bw_lorenzen_spectrum')

if __name__=='__main__':main()
