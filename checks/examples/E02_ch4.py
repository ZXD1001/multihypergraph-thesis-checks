"""Three-qutrit joint effects, verification spectrum, and sample-count bounds."""

from pathlib import Path
import os
from itertools import product
import numpy as np
import hashlib, json, math

source_root=os.environ.get('THESIS_AUDIT_SOURCE_ROOT')
root=Path(source_root) if source_root else None
names=['data/chap04.tex','data/appendix_w2.tex','data/cnabstract.tex','data/enabstract.tex','data/chap01.tex','data/chap07.tex']
before={n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in names} if root is not None else {}

def kron(vs):
    out=np.array([1.0+0j])
    for v in vs: out=np.kron(out,v)
    return out

d=3
labels=list(product(range(d),repeat=3))
w=np.exp(2j*np.pi/d)
F=np.array([[w**(-x*j)/np.sqrt(d) for j in range(d)] for x in range(d)])
D=np.diag([w**(x[0]*x[1]+2*x[0]*x[1]**2) for x in labels])
U=D@np.kron(np.kron(F,F),F)
Ps=[U@np.diag([int(j[v]==0) for j in labels])@U.conj().T for v in range(3)]
Es=[]; residuals=[]
for A in [(0,2),(1,)]:
    B=tuple(v for v in range(3) if v not in A)
    # Check every second-layer outcome, not only accepting outcome.
    for ja in product(range(d),repeat=len(A)):
        E=np.zeros((27,27),complex)
        for yy in product(range(d),repeat=len(B)):
            y=dict(zip(B,yy)); j=dict(zip(A,ja)); vecs=[]
            for v in range(3):
                if v in B: vecs.append(np.eye(d)[:,y[v]])
                else:
                    phase=lambda x: (y[1]+2*y[1]**2)*x if v==0 else (y[0]*(x+2*x*x) if v==1 else 0)
                    vecs.append(np.array([w**(phase(x)-j[v]*x)/np.sqrt(d) for x in range(d)]))
            vv=kron(vecs); E+=np.outer(vv,vv.conj())
        expected=U@np.diag([int(all(label[v]==j[v] for v in A)) for label in labels])@U.conj().T
        residuals.append(float(np.linalg.norm(E-expected,ord=2)))
        if all(x==0 for x in ja): Es.append(E)
eig=np.linalg.eigvalsh((Es[0]+Es[1])/2)
unique,count=np.unique(np.round(eig,10),return_counts=True)
eps=.01; delta=.05; nu=.5
exact=math.ceil(math.log(delta)/math.log1p(-eps*nu)); simple=math.ceil(math.log(1/delta)/(eps*nu))
assert exact==598 and simple==600
assert max(residuals)<1e-12
assert list(zip(unique,count))==[(0.,16),(.5,10),(1.,1)]
assert math.exp(exact*math.log1p(-eps*nu))<=delta
assert math.exp((exact-1)*math.log1p(-eps*nu))>delta

after={n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in names} if root is not None else {}
assert before==after
result={'checked':'3-qutrit actual joint effects for all second-layer outcomes; weighted spectrum; sample ceiling',
'max_operator_norm_residual':max(residuals),'equal_weight_spectrum':{str(v):int(n) for v,n in zip(unique,count)},
'exact_N':exact,'simplified_N':simple,'allpass_bound_N598':math.exp(exact*math.log1p(-eps*nu)),
'allpass_bound_N597':math.exp((exact-1)*math.log1p(-eps*nu)),
'correlated_counterexample_allpass_lower_bound':1-eps,'source_hashes_unchanged':True,'source_sha256':before}
print(json.dumps(result,ensure_ascii=False,indent=2))
