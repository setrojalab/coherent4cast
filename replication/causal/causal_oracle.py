"""Independent dense matrix-exponential g-formula, not recurcast uniformization."""
from pathlib import Path
import csv
import numpy as np
from scipy.linalg import expm
out=Path(__file__).resolve().parent/'results'
out.mkdir(parents=True,exist_ok=True)
states=[(0,0),(0,1),(1,0),(1,1)]
results=[]
for w in [0,1]:
    values=[]
    for x in [0,1]:
        Q=np.zeros((5,5)); reward=np.zeros((5,4)); terminal=np.zeros((5,4));terminal[4,2]=1
        for s,(a,b) in enumerate(states):
            rates=np.array([.8,.35,.4])*np.exp([.25*x-.35*w+.15*a+.10*b,
                .15*x+.10*w+.10*a+.20*b,.30*x-.25*w+.25*a+.20*b])
            nexts=[states.index((min(1,a+1),b)),states.index((a,min(1,b+1))),4]
            for k,t in enumerate(nexts):Q[s,t]+=rates[k];Q[s,s]-=rates[k]
            reward[s,:2]=rates[:2];reward[s,3]=1
        augmented=np.block([[Q,reward],[np.zeros((4,5)),np.zeros((4,4))]])
        transition=expm(augmented)
        v=transition[:5,:5]@terminal+transition[:5,5:]
        assert np.allclose(Q.sum(axis=1),0) and np.allclose(transition[:5,:5].sum(axis=1),1)
        values.append(v[0])
    mean=np.mean(values,axis=0)
    for target,value in zip(['burden_1','burden_2','risk','rmst'],mean):
        results.append(dict(arm=w,target=target,oracle=value))
with (out/'dense_gformula_oracle.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['arm','target','oracle']);writer.writeheader();writer.writerows(results)
print(results)
