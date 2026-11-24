"""Finite examples separating branch, channel, and classical privacy properties."""

import numpy as np,json
q=5
X=np.roll(np.eye(q),1,axis=0)
phi=np.eye(q).reshape(-1)/np.sqrt(q)
a=np.outer(phi,phi.conj())
U=np.kron(np.eye(q),X)
J=(a+U@a@U.conj().T)/2
w=np.exp(2j*np.pi/q)
M=np.array([w**(z**3) for z in range(q)])/np.sqrt(q)
T=np.diag([w**(z**3) for z in range(q)])
V=[T@np.linalg.matrix_power(X,t)@T.conj().T for t in range(q)]
mu=np.outer(M,M.conj())
F=np.array([[w**(z*k) for k in range(q)] for z in range(q)])/np.sqrt(q)
e=T@F
O=np.outer(e[:,0],e[:,1].conj())
twirl=lambda x:sum(u@x@u.conj().T for u in V)/q
key_marginals=[]
cipher_marginals=[]
recovery=[]
for s in range(q):
    p=np.zeros((q,q))
    for r in range(q):p[r,(s-r)%q]=1/q
    key_marginals.append(p.sum(axis=1));cipher_marginals.append(p.sum(axis=0))
    recovery.append(sum(p[r,t] for r in range(q) for t in range(q) if (r+t)%q==s))
d={'noncanonical_choi_rank':int(np.linalg.matrix_rank(J,tol=1e-10)),'nonzero_choi_eigenvalues':np.linalg.eigvalsh(J)[-2:].tolist(),'choi_reference_marginal_error':float(np.max(np.abs(J.reshape(q,q,q,q).trace(axis1=1,axis2=3)-np.eye(q)/q))),'unrecorded_twirl_magic_fixed_error':float(np.max(np.abs(twirl(mu)-mu))),'unrecorded_twirl_offdiagonal_norm':float(np.linalg.norm(twirl(O))),'key_view_max_difference':float(np.max(np.ptp(np.array(key_marginals),axis=0))),'cipher_view_max_difference':float(np.max(np.ptp(np.array(cipher_marginals),axis=0))),'all_secret_recovery_probabilities':recovery}
print(json.dumps(d,indent=2))
assert d['noncanonical_choi_rank']==2
assert all(abs(p-1)<1e-10 for p in recovery)
assert d['unrecorded_twirl_offdiagonal_norm']<1e-10
