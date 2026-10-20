import pickle, itertools, sys, time
import numpy as np
from w9b_k3search import chain_inst
from w9b_pend import *
good=pickle.load(open('/tmp/w9b_k3good.pkl','rb'))
import re,ast
fam_by_idx={}
for line in open('/tmp/k3fam.log'):
    m=re.match(r'(\d+) sig .*family \(member,list\): (.*)$',line.strip())
    if m: fam_by_idx[int(m.group(1))]=ast.literal_eval(m.group(2))
def lift_family(fam4,N):
    """pattern from N=4 family: member index m=3*j+t. 2-lists at type t for all j; singleton lists at blocks 0/N-1 ends.  Return mapping function to build N-family."""
    types=sorted({(m%3,tuple(L)) for m,L in fam4 if len(L)==2})
    ts=[t for t,_ in types]; assert len(ts)==1, types
    t2=ts[0]
    out=[(3*j+t2,(0,1)) for j in range(N)]
    for m,L in fam4:
        if len(L)==1:
            j,t=divmod(m,3)
            # block index: first end (j<2) -> block 0 ; last end -> block N-1
            jj=0 if j<2 else N-1
            out.append((3*jj+t,tuple(L)))
    return out
res=[]
N=int(sys.argv[1])
t0=time.time()
for idx,(R,sig,cfg,outm,inm,eps,ns,k2) in enumerate(good):
    if idx not in fam_by_idx: continue
    fam=lift_family(fam_by_idx[idx],N)
    I=chain_inst(N,cfg,outm,inm,eps)
    S=enum_sigma(I,limit=10**6,drop_far=True)
    Sg=np.array(S,dtype=np.int8)
    # verify family inconsistent & minimal on Sigma
    def feas(fm):
        ok=np.ones(len(Sg),bool)
        for m,L in fm: ok&=np.isin(Sg[:,m],list(L))
        return ok.any()
    incons=not feas(fam); minimal=all(feas([f for f in fam if f!=g]) for g in fam)
    # pendant choices for the apart-type guards
    best=None
    nz=[k for k,(m,L) in enumerate(fam) if LIST2ROLE[tuple(sorted(L))][1]==0]
    for ecs in itertools.product((0,1),repeat=min(len(nz),4)):
        ec=[0]*len(fam)
        # use same choice pattern for all 2-list guards: ecs[0]; singletons choose by ecs[1..]
        for k in nz: ec[k]=ecs[0] if len(fam[k][1])==2 else ecs[min(1,len(ecs)-1)]
        I2,T=add_pendants(I,fam,ec)
        if not I2.valid(): continue
        r=min_certificate_size(I2,limit_budget=12)
        mc=r[0] if r else None
        if best is None or (mc or 0)>(best[0] or 0): best=(mc,tuple(ec),len(T))
    res.append((idx,sig,outm,inm,incons,minimal,best))
    print(idx,'sig',sig,'out/in',outm,inm,'family inconsistent/minimal on Sigma:',incons,minimal,'best min-certificate over pendant end choices:',best,'intended |T|',best[2] if best else None,round(time.time()-t0),flush=True)
