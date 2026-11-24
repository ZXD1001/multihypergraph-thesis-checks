"""Finite-field certificates, direction ranks, and gate-sequence checks."""

from itertools import product
import json,numpy as np
q=5
rows=list(product(range(q),repeat=3))
B=np.zeros((125,25),dtype=int)
J=np.zeros((125,5),dtype=int)
y=[];g={}
for r,(a,b,c) in enumerate(rows):
    vals=[(a+b*i+c*i*i)%q for i in range(q)]
    for i,v in enumerate(vals):B[r,i*q+v]=1
    J[r,c]=1
    y.append(((1 if a==0 else 0)*[1,3,1,0,0][b]*(4-c))%q)
    assert sum(vals)%q==0
    g[a,b,c]=sum(vals)//q
h={}
for b,c in product(range(q),repeat=2):
    h[0,b,c]=0
    for a in range(4):h[a+1,b,c]=(h[a,b,c]+c**3-g[a,b,c])%q
    assert (h[0,b,c]-h[4,b,c]-c**3+g[4,b,c])%q==0
v=np.array([-(h[(a+1)%q,b,c]-h[a,b,c]) for a,b,c in rows])
y=np.array(y)
res={'rows':len(rows),'B_shape':list(B.shape),'B_nonzeros':int(np.count_nonzero(B)),
'D23_cyclic_checks':25,'D27_yB_zero_mod5':bool(np.all((y@B)%q==0)),
'D27_yJ_zero_mod5':bool(np.all((y@J)%q==0)),'D27_y_minus_delta_h_mod5':int((y@v)%q),
'D27_yB_exact_integer_zero':bool(np.all(y@B==0)),'D27_yB_integer_max':int(np.max(y@B))}
assert res['D27_yB_zero_mod5'] and res['D27_yJ_zero_mod5'] and res['D27_y_minus_delta_h_mod5']==1 and not res['D27_yB_exact_integer_zero']
def rank5(A):
    A=np.array(A,dtype=int)%q; r=0
    for c in range(A.shape[1]):
        piv=next((i for i in range(r,len(A)) if A[i,c]),None)
        if piv is None:continue
        A[[r,piv]]=A[[piv,r]];A[r]=(A[r]*pow(int(A[r,c]),-1,q))%q
        for i in range(len(A)):
            if i!=r:A[i]=(A[i]-A[i,c]*A[r])%q
        r+=1
        if r==len(A):break
    return r
res['polynomial_direction_ranks']=[rank5([[c**d for d in range(1,k+1)] for c in range(q)]) for k in (2,3,4)]
res['logical_counts_from_direction_ranks']=[q**r for r in res['polynomial_direction_ranks']]
light=[((b+t)*c**3-b*c**3-t*c**3)%q for a,b,c,s,t in product(range(q),repeat=5)]
res['light_shift_identity_checks']=len(light);res['light_shift_identity_all_zero']=all(x==0 for x in light)
def matops(n,ops):
    A=np.eye(n,dtype=int)
    for op in ops:
        if op[0]=='add':_,dst,src,w=op;A[dst]=(A[dst]+w*A[src])%q
        elif op[0]=='scale':_,i,w=op;A[i]=(w*A[i])%q
        else:_,i,j=op;A[[i,j]]=A[[j,i]]
    return A
E=matops(5,[('add',3,0,1),('add',3,1,3),('add',3,2,4),('add',4,0,1),('add',4,1,4),('add',4,2,1),('add',1,2,1),('add',1,0,1),('scale',2,2),('add',2,1,2),('add',2,0,-1)])
U=matops(3,[('add',1,2,1),('scale',2,4),('add',2,1,4),('scale',1,2),('add',2,0,3),('add',1,0,1)])
T=matops(3,[('add',1,2,2),('add',0,2,3),('add',0,1,3),('swap',1,2),('add',2,0,1),('swap',0,1)])
res['D10_E5_match']=bool(np.array_equal(E,[[1,0,0,0,0],[1,1,1,0,0],[1,2,4,0,0],[1,3,4,1,0],[1,4,1,0,1]]))
res['D11_U012_match']=bool(np.array_equal(U,[[1,0,0],[1,2,2],[3,4,3]]))
res['D12_T34_match']=bool(np.array_equal(T,[[0,0,1],[1,3,4],[1,4,1]]))
assert all(res[k] for k in ['D10_E5_match','D11_U012_match','D12_T34_match','light_shift_identity_all_zero'])
print(json.dumps(res,indent=2))
