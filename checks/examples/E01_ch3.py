"""Finite enumeration of phase-state families and composite-dimension null polynomials."""

from itertools import product
from math import prod
import json

def image_set(d,n,predicate):
    inputs=list(product(range(d),repeat=n))
    exponents=[a for a in product(range(d),repeat=n) if any(a) and predicate(a)]
    out={tuple(0 for _ in inputs)}
    for a in exponents:
        g=tuple(prod(pow(xi,ai,d) for xi,ai in zip(x,a))%d for x in inputs)
        out={tuple((vi+c*gi)%d for vi,gi in zip(v,g)) for v in out for c in range(d)}
    return out,len(exponents)

results=[]
for d,n in [(2,2),(3,2),(4,1),(4,2),(6,1)]:
    families={}
    for name,pred in [('G',lambda a:sum(x!=0 for x in a)==2 and all(x in (0,1) for x in a)),('MG',lambda a:sum(x!=0 for x in a)==2),('HG',lambda a:all(x in (0,1) for x in a)),('MHG',lambda a:True)]:
        im,k=image_set(d,n,pred)
        families[name]=im
        assert d**k%len(im)==0
    assert families['G']==families['MG']&families['HG']
    assert families['G']<=families['MG']<=families['MHG']
    assert families['HG']<=families['MHG']
    results.append({'d':d,'n':n,'distinct_states':{k:len(v) for k,v in families.items()},'intersection_check':True})

composites=[]
for d in range(4,61):
    p=next((p for p in range(2,d) if d%p==0),None)
    if p is None: continue
    assert p<d
    coeff=[1]
    for j in range(p):
        nxt=[0]*(len(coeff)+1)
        for i,c in enumerate(coeff):
            nxt[i]-=j*c
            nxt[i+1]+=c
        coeff=nxt
    coeff=[(d//p*c)%d for c in coeff]
    assert coeff[0]==0 and coeff[-1]==d//p and any(coeff)
    assert all(sum(c*pow(x,i,d) for i,c in enumerate(coeff))%d==0 for x in range(d))
    composites.append(d)
print(json.dumps({'family_tests':results,'explicit_null_polynomial_composites_d_le_60':composites},ensure_ascii=False,indent=2))
