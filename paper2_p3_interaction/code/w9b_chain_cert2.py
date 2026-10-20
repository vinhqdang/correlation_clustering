import sys, time
from w9b_chain_cert import *
for N in (6,8,10):
    I,T=build(N); far=set(I.far_pairs())
    t0=time.time()
    c_full,w=leaf_cost(I,T); c_empty,_=leaf_cost(I,[])
    bad=[e for e in T if leaf_cost(I,[p for p in T if p!=e])[0]!=len(I.H)]
    o_full=opt_csp_cpsat(I,T)
    print('N',N,'members',len(I.H),'n',I.n,'|T|',len(T),'F(T)-|H|=',c_full-len(I.H),'F(empty)-|H|=',c_empty-len(I.H),'unnecessary elements of T:',bad,'cpsat F(T)-|H|=',o_full-len(I.H),'elim width',w,round(time.time()-t0,1),flush=True)
