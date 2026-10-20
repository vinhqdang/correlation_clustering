import pickle, numpy as np, time, sys
from w9b_k3search import *
good=pickle.load(open('/tmp/w9b_k3good.pkl','rb'))
out=[]
for idx,(R,sig,cfg,outm,inm,eps,ns,k2) in enumerate(good):
    I=chain_inst(4,cfg,outm,inm,eps)
    S=enum_sigma(I,limit=10**5); Sg=np.array(S,dtype=np.int8)
    k,fam=max_core(Sg,maxk=9,return_family=True)
    out.append((idx,k,fam))
    print(idx,'sig',sig,'out/in',outm,inm,'|S4|',len(S),'kappa2(N=4)=',k,'family (member,list):',sorted(fam),flush=True)
pickle.dump(out,open('/tmp/w9b_k3fam.pkl','wb'))
