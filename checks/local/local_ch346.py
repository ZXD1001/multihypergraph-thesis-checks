from pathlib import Path
import os
from fractions import Fraction
from math import lcm
from itertools import product
import json, hashlib, math
import numpy as np

OUT=(Path(os.environ['THESIS_WEB_AUDIT_RESULTS'])/'local_ch346.json')
q=5
rows=list(product(range(q),repeat=3))
targets=list(product(range(q),repeat=4))
B=np.zeros((125,25),dtype=np.int64)
C=np.zeros((125,15),dtype=np.int64)
for z,(a,b,c) in enumerate(rows):
    for i in range(q):
        v=(a+b*i+c*i*i)%q
        B[z,5*i+v]=1
        C[z,3*i:3*i+3]=[1,v,v*v%q]
P=np.array([[sum(t[j-1]*pow(c,j,q) for j in range(1,5))%q for t in targets] for a,b,c in rows],dtype=np.int64)

def rational_kernel(A):
    R=[[Fraction(int(x)) for x in row] for row in A]
    piv=[];r=0
    for j in range(len(R[0])):
        k=next((k for k in range(r,len(R)) if R[k][j]),None)
        if k is None:continue
        R[r],R[k]=R[k],R[r];s=R[r][j];R[r]=[x/s for x in R[r]]
        for k in range(len(R)):
            if k!=r and R[k][j]:
                s=R[k][j];R[k]=[x-s*y for x,y in zip(R[k],R[r])]
        piv.append(j);r+=1
        if r==len(R):break
    N=[]
    for j in range(len(R[0])):
        if j in piv:continue
        v=[Fraction(0)]*len(R[0]);v[j]=1
        for i,k in enumerate(piv):v[k]=-R[i][j]
        d=lcm(*(x.denominator for x in v))
        N.append([int(d*x) for x in v])
    return np.array(N,dtype=np.int64),piv

def modular_elimination(A):
    R=A.copy()%q;E=np.eye(len(A),dtype=np.int64);piv=[];r=0
    for j in range(R.shape[1]):
        candidates=np.flatnonzero(R[r:,j])
        if not len(candidates):continue
        k=r+int(candidates[0]);R[[r,k]]=R[[k,r]];E[[r,k]]=E[[k,r]]
        s=pow(int(R[r,j]),-1,q);R[r]=(R[r]*s)%q;E[r]=(E[r]*s)%q
        for k in range(len(R)):
            if k!=r and R[k,j]:
                s=int(R[k,j]);R[k]=(R[k]-s*R[r])%q;E[k]=(E[k]-s*E[r])%q
        piv.append(j);r+=1
    assert np.array_equal((E@A)%q,R)
    return R,E,piv

N,rp=rational_kernel(B.T)
assert np.all(N@B==0)
R,E,piv=modular_elimination(B)
RC,EC,pc=modular_elimination(C)
NP=(N@P)%q;EP=(E@P)%q;ECP=(EC@P)%q
functions={'0':lambda a,b,c:0,'bc3':lambda a,b,c:b*c**3,'ac3+bc4':lambda a,b,c:a*c**3+b*c**4,'ac3+a4':lambda a,b,c:a*c**3+a**4}
stats={};cert75=[];positive_total=0;negative_total=0
for name,h in functions.items():
    lu=set();lc=set();per=[]
    for s,t in product(range(q),repeat=2):
        delta=np.array([(h((a+s)%q,(b+t)%q,c)-h(a,b,c))%q for a,b,c in rows],dtype=np.int64)
        residual=(EP-(E@delta)[:,None])%q
        good=np.all(residual[len(piv):]==0,axis=0)
        neg=(NP-(N@delta)[:,None])%q
        rejected=np.any(neg!=0,axis=0)
        assert np.all(good|rejected), (name,s,t,'unresolved')
        assert not np.any(good&rejected)
        indices=np.flatnonzero(good)
        U=np.zeros((25,len(indices)),dtype=np.int64)
        for r,j in enumerate(piv):U[j]=residual[r,indices]
        assert np.all((B@U-P[:,indices]+delta[:,None])%q==0)
        cl=np.all(((ECP-(EC@delta)[:,None])%q)[len(pc):]==0,axis=0)
        lu.update(map(int,indices));lc.update(map(int,np.flatnonzero(cl)))
        positive_total+=int(good.sum());negative_total+=int(rejected.sum())
        per.append({'shift':[s,t],'lu_targets':int(good.sum()),'lc_targets':int(cl.sum())})
        if name!='ac3+bc4':
            target=(0,0,1,0) if name=='0' else (0,0,0,1)
            ix=targets.index(target);wi=int(np.flatnonzero(neg[:,ix])[0]);y=N[wi]
            assert np.all(y@B==0) and int(y@(P[:,ix]-delta))%q!=0
            cert75.append({'h':name,'target_coefficients_c1_to_c4':list(target),'shift':[s,t],'integer_kernel_row':wi,'residual_mod5':int(y@(P[:,ix]-delta))%q})
    stats[name]={'lu_count':len(lu),'lc_count':len(lc),'per_shift':per}
assert [(stats[h]['lc_count'],stats[h]['lu_count']) for h in functions]==[(25,25),(125,125),(625,625),(25,125)]
assert len(cert75)==75

# Independent full phase back-substitution for all 3 original generators, all 125 rows.
algorithm_A=[]
for target,shift,gamma2 in [((0,0,1,0),(1,0),0),((0,1,2,0),(2,0),1),((0,0,0,1),(0,1),0)]:
    s,t=shift; count=0
    for a,b,c in rows:
        ap=(a-s)%q;bp=(b-t)%q
        physical=sum(-gamma2*((ap+bp*i+c*i*i)%q)**2 for i in range(q))
        phase=functions['ac3+bc4'](a,b,c)-functions['ac3+bc4'](ap,bp,c)+physical
        assert (phase-sum(target[j-1]*c**j for j in range(1,5)))%q==0
        count+=1
    algorithm_A.append({'target_coefficients_c1_to_c4':list(target),'shift':list(shift),'quadratic_coefficient':gamma2,'exact_row_checks':count})

# Sample-count thresholds are checked with integer powers, not floating logarithms.
samples=[]
for nu,Ns in [(Fraction(1,10),2995),(Fraction(1,4),1197),(Fraction(1,2),598),(Fraction(1,3),898)]:
    x=1-nu/Fraction(100);num,den=x.numerator,x.denominator
    assert 20*pow(num,Ns)<=pow(den,Ns)
    assert 20*pow(num,Ns-1)>pow(den,Ns-1)
    samples.append({'nu':str(nu),'N':Ns,'integer_N_and_previous_verified':True})

result={'origin':'exact finite checks of phase libraries, compilation, and sample counts','scope':'q=5,k=3,m=1, ordered support 0..4, four stated h, all 25 shifts and 625 phase targets; continuous-LU exclusions certified by exact integer kernel rows and feasible phases by mod-5 positive certificates; not a general-parameter proof',
 'B_shape':list(B.shape),'B_integer_rank':len(rp),'B_mod5_rank':len(piv),'C_mod5_rank':len(pc),'integer_kernel_shape':list(N.shape),'integer_kernel_max_abs':int(np.max(np.abs(N))), 'integer_kernel_saturation':'not claimed or needed for these results: every examined case independently has either a checked positive or a checked negative certificate',
 'statistics':stats,'positive_shift_target_cases':positive_total,'negative_shift_target_cases':negative_total,'certificate_75':cert75,'integer_kernel_rows':N.tolist(),'algorithm_A':algorithm_A,'sample_counts':samples}
OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps({'result':str(OUT),'counts':{h:{k:v for k,v in x.items() if k!='per_shift'} for h,x in stats.items()},'positive_cases':positive_total,'negative_cases':negative_total,'negative_certificates':len(cert75),'algorithm_A_row_checks':sum(x['exact_row_checks'] for x in algorithm_A),'B_integer_rank':len(rp),'B_mod5_rank':len(piv),'integer_kernel_max_abs':int(np.max(np.abs(N)))},ensure_ascii=False,indent=2))
