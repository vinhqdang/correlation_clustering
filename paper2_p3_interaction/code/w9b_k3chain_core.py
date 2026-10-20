import numpy as np, itertools, sys, time
from w9b_k3chain import *
def analyse(N):
    I=chain(N,(0,0,0,0))
    S=enum_sigma(I,limit=10**6)
    Sg=np.array(S,dtype=np.int8)
    m=len(I.H)
    # members order: for j: A_j=3j, B_j=3j+1, C_j=3j+2
    def sat(mask_members):
        pass
    rows=Sg
    def feasible(lists):   # lists: dict member -> allowed set
        ok=np.ones(len(rows),bool)
        for i,L in lists.items(): ok&=np.isin(rows[:,i],list(L))
        return ok.any()
    res=[]
    for a in (0,1,2):
        for b in (0,1,2):
            fam={3*0:{a},3*(N-1)+1:{b}}
            for j in range(N): fam[3*j+2]={0,1}
            if not feasible(fam):
                # minimality
                minimal=all(feasible({k:v for k,v in fam.items() if k!=g}) for g in fam)
                res.append((a,b,minimal))
    return I,len(S),res
for N in range(1,8):
    t0=time.time()
    I,s,res=analyse(N)
    print('N',N,'members',len(I.H),'|Sigma|',s,'inconsistent endpoint choices (a,b,minimal):',res,round(time.time()-t0,1),flush=True)
