"""Choi and discrete-Wigner diagnostics for a two-branch five-dimensional channel."""

import json,numpy as np
q=5;w=np.exp(2j*np.pi/q)
I=np.eye(q);X=np.roll(I,1,axis=0)
phi=np.eye(q).reshape(-1)/np.sqrt(q)
phi1=np.kron(I,X)@phi
J=(np.outer(phi,phi.conj())+np.outer(phi1,phi1.conj()))/2
M=np.array([w**(z**3) for z in range(q)])/np.sqrt(q)
rho=np.outer(M,M.conj());out=(rho+X@rho@X.conj().T)/2
inv2=3
W=np.array([[(sum(w**(p*y)*out[(x-inv2*y)%q,(x+inv2*y)%q] for y in range(q))/q).real for p in range(q)] for x in range(q)])
res={'q':q,'ancilla_probabilities':[.5,.5,0,0,0],'choi_eigenvalues':np.linalg.eigvalsh(J).round(14).tolist(),'choi_rank':int(np.linalg.matrix_rank(J,tol=1e-12)),'choi_trace':float(np.trace(J).real),'output_wigner_min':float(W.min()),'output_wigner_sum':float(W.sum()),'max_branch_effect_error':max(float(np.max(abs(K.conj().T@K-I/2))) for K in [I/np.sqrt(2),X/np.sqrt(2)])}
print(json.dumps(res,indent=2))
