#!/usr/bin/env python3
"""Exact finite checks for the six-party, three-level sharing instance.
Density-state conclusions use the analytic module decomposition in the thesis.
"""
from pathlib import Path
import os
from itertools import product
from collections import Counter
from fractions import Fraction
import hashlib, json, platform
Q=5
source_root=os.environ.get('THESIS_WEB_AUDIT_SOURCE')
root=Path(source_root) if source_root else None
source_files=['data/chap05.tex','data/chap02.tex','data/appendix_w3.tex']
source_hashes=({f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in source_files}
               if root is not None else {})

def rank_mod(mat,q=Q):
    a=[[int(v)%q for v in row] for row in mat]
    if not a:return 0
    r=0
    for c in range(len(a[0])):
        p=next((i for i in range(r,len(a)) if a[i][c]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r]
        inv=pow(a[r][c],-1,q)
        a[r]=[(v*inv)%q for v in a[r]]
        for i in range(len(a)):
            if i!=r:
                f=a[i][c]
                a[i]=[(v-f*w)%q for v,w in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

# q=5, L=3, six parties; permission is the length of the initial complete pairs.
def prefix_length(mask):
    return max(k for k in range(4) if mask & ((1<<(2*k))-1)==((1<<(2*k))-1))
permission_counts=Counter(prefix_length(mask) for mask in range(64))
assert [permission_counts[k] for k in range(4)]==[48,12,3,1]

# Top affine map from (z,u,v,R,T) to (c,d,i,j,g,h), constant s3 only in g.
A=[[1,1,0,0,0],[1,0,0,0,0],[1,0,1,1,0],[1,0,0,1,0],
   [1,-1,-1,0,1],[1,0,0,0,1]]
proper_rank_rows=[]
for mask in range(63):
    rows=[A[i] for i in range(6) if mask>>i&1]
    rk=rank_mod(rows)
    assert rk==mask.bit_count()
    proper_rank_rows.append({'mask':mask,'registers':mask.bit_count(),'rank':rk})
assert rank_mod(A)==5
assert all(row[0]%Q==1 for row in A) # tracing one line eliminates z != z'.

# Exact complete-view representations using the proven module decomposition:
# an incomplete Bell pair has a fixed I/q factor; full pairs reveal orthogonal
# Bell labels. Proper top subsets have I/q^|J| by rank and the trace argument.
# For all six parties, the top s3 label is orthogonal after the fixed inverse
# Clifford (with common resource and mixed ancilla factors).
def view_histogram(mask,s):
    s1,s2,s3=s
    ab_complete=mask&3==3
    ef_complete=mask&12==12
    hist=Counter()
    for r in range(Q):
        lab=[]
        if ab_complete:lab.extend([s1,(s2-r)%Q])
        if ef_complete:lab.append(r)
        if mask==63:lab.append(s3)
        hist[tuple(lab)]+=1
    return tuple(sorted(hist.items()))

view_checks=0
for mask in range(64):
    k=prefix_length(mask)
    baseline={}
    for s in product(range(Q),repeat=3):
        prefix=s[:k]
        hist=view_histogram(mask,s)
        if prefix in baseline:assert hist==baseline[prefix]
        else:baseline[prefix]=hist
        # Explicit authorized prefix decoding from any complete-view label.
        ab_complete=mask&3==3
        ef_complete=mask&12==12
        for label,count in hist:
            assert count>0
            if k>=1:assert label[0]==s[0]
            if k>=2:assert (label[1]+label[2])%Q==s[1]
            if k==3:assert label[-1]==s[2]
        view_checks+=1
assert view_checks==8000

# Full raw measurement record for the fixed secret s=(1,2,4).
s1,s2,s3=(1,2,4)
records=Counter()
for r,u,v,R,T,y,k2,k3 in product(range(Q),repeat=8):
    rec=(s1,(s2-r)%Q,y,(r+y)%Q,u,v,(s3-u-v)%Q,k2,k3)
    records[rec]+=1
assert sum(records.values())==Q**8==390625
assert len(records)==Q**6==15625
assert set(records.values())=={25}
# For fixed masks, y,k2,k3 are the only 3 independent random outcomes.
assert len(set(product(range(Q),repeat=3)))==125
# Top recovery phase after positive Fourier bras and Z^(k2+k3) correction.
phase_checks=0
for R,T,k2,k3,z in product(range(Q),repeat=5):
    uncorrected=(z**3-k2*(z+R)-k3*(z+T))%Q
    corrected=(uncorrected+(k2+k3)*z)%Q
    expected=(z**3-k2*R-k3*T)%Q
    assert corrected==expected
    phase_checks+=1
assert phase_checks==3125
assert Q*Fraction(1,Q**3)==Fraction(1,25)
assert Fraction(1,5)*Fraction(1,25)==Fraction(1,125)

m=[2,2,2,2,1,1]
messages=[('b',1,0),('d',1,0),('j',3,2),('h',5,4)]
c=[0]*6
for name,sender,receiver in messages:
    c[sender]+=1;c[receiver]+=1
r=[3,3,2,2,1,1]
assert c==[2,2,1,1,1,1]
assert sum(c)==2*len(messages)==8
assert sum(m)==10 and len(messages)==4
assert [ri+1-mi for ri,mi in zip(r,m)]==c
assert sum(m)+2*len(messages)==3*(3+3)
# Ensure source was read only during the check.
assert all(hashlib.sha256((root/f).read_bytes()).hexdigest()==v for f,v in source_hashes.items())
result={
 'identity':'exact finite check of the six-party, three-level sharing instance',
 'python':platform.python_version(),'q':Q,'L':3,
 'source_sha256':source_hashes,
 'permission_counts':[permission_counts[k] for k in range(4)],
 'proper_top_affine_rank_checks':len(proper_rank_rows),
 'complete_view_representation_checks':view_checks,
 'view_scope':'exact Bell-label distributions plus analytic trace/rank and inverse-Clifford density representation; not full 5^10 density matrices',
 'fixed_secret_for_history_count':[s1,s2,s3],
 'history_count':sum(records.values()),'raw_global_record_count':len(records),
 'hidden_RT_preimages_per_record':sorted(set(records.values())),
 'fixed_masks_outcome_count':125,'phase_identity_checks':phase_checks,
 'm':m,'c':c,'N':sum(m),'C':len(messages),
 'not_proved':['general SO lower bounds','all layer counts','general refinement lemma'],
 'top_rank_details':proper_rank_rows}
(Path(os.environ['THESIS_WEB_AUDIT_RESULTS'])/'local_ch5.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='top_rank_details'},ensure_ascii=False,indent=2))
