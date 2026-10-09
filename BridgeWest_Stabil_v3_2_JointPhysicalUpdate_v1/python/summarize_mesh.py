from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.io import loadmat,savemat
r=Path(__file__).resolve().parents[1]
A=loadmat(r/'results/python_model.mat');B=loadmat(r/'results/refined_model.mat')
lookup={tuple(p.round(7)):i for i,p in enumerate(B['nodes'])}
pairs=[(i,lookup[tuple(p.round(7))]) for i,p in enumerate(A['nodes']) if tuple(p.round(7)) in lookup]
ai,bi=np.array(pairs).T;ia=(6*ai[:,None]+np.arange(3)).ravel();ib=(6*bi[:,None]+np.arange(3)).ravel()
w=np.maximum(A['M'].diagonal()[ia],0);U=A['modes'][ia]*np.sqrt(w)[:,None];V=B['modes'][ib]*np.sqrt(w)[:,None]
mac=(U.T@V)**2/(np.sum(U*U,axis=0)[:,None]*np.sum(V*V,axis=0)[None,:]);fa=A['frequencies_Hz'][0];fb=B['frequencies_Hz'][0]
rows=[]
for j,label in [(2,'first bending'),(7,'first torsion'),(108,'second bending candidate')]:
 k=int(np.argmax(mac[j]));change=100*(fb[k]/fa[j]-1)
 rows.append(dict(family=label,base_fe_index=j+1,refined_fe_index=k+1,base_frequency_Hz=fa[j],refined_frequency_Hz=fb[k],frequency_change_percent=change,common_node_translation_MAC=mac[j,k],common_nodes=len(pairs),status='stable in this refinement' if abs(change)<2 and mac[j,k]>.9 else 'not converged / assignment unstable'))
with (r/'results/mesh_tracking.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
# Synchronise corrected provenance text; numerical inputs/results unchanged.
c=json.loads((r/'config/model.json').read_text())
for name in ['python_model.mat','refined_model.mat']:
 d=loadmat(r/'results'/name);cc=json.loads(str(d['config_json'][0]));cc['geometry_status']=c['geometry_status'];d['config_json']=json.dumps(cc)
 savemat(r/'results'/name,{k:v for k,v in d.items() if not k.startswith('__')},do_compression=True)
# Retain the unmodified source identities for all reused routines.
(r/'evidence/stabil_source_hashes.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (r/'dependencies/Stabil').glob('*.m')},indent=2))
print(json.dumps(rows,indent=2))
