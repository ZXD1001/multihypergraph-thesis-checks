"""Finite cross-operator privacy checks for two polynomial-code instances."""

import itertools, json
import numpy as np

results=[]
for q,k in [(5,2),(5,3)]:
    n=2*k-1
    rng=np.random.default_rng(20261025+k)
    # A general phase table, rather than only a special low-degree phase.
    phase=rng.integers(0,q,size=(q,)*k)
    code=np.zeros((q**n,q),dtype=complex)
    omega=np.exp(2j*np.pi/q)
    for c in range(q):
        for a in itertools.product(range(q),repeat=k-1):
            coef=a+(c,)
            ys=tuple(sum(coef[j]*(x**j) for j in range(k))%q for x in range(n))
            code[np.ravel_multi_index(ys,(q,)*n),c]=omega**int(phase[coef])/np.sqrt(q**(k-1))
    orth_error=float(np.max(np.abs(code.conj().T@code-np.eye(q))))
    max_error=0.0
    cases=0
    tensors=[code[:,c].reshape((q,)*n) for c in range(q)]
    for size in range(k):
        for A in itertools.combinations(range(n),size):
            B=tuple(i for i in range(n) if i not in A)
            mats=[v.transpose(A+B).reshape(q**len(A),q**len(B)) for v in tensors]
            for c,cp in itertools.product(range(q),repeat=2):
                actual=mats[c]@mats[cp].conj().T
                expected=(np.eye(q**size)/q**size) if c==cp else np.zeros_like(actual)
                max_error=max(max_error,float(np.max(np.abs(actual-expected))))
                cases+=1
    assert orth_error<1e-12 and max_error<1e-12
    results.append({'q':q,'k':k,'n':n,'m':1,'phase':'deterministic random full table','seed':20261025+k,'cross_operator_checks':cases,'max_entry_error':max_error,'isometry_max_entry_error':orth_error})
print(json.dumps(results,indent=2))
