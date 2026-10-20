"""partition-variant (wave 8 definition: far pairs outside T free, transitivity kept) check of the K3-chain certificate."""
import time, sys
from w9b_chain_cert import *
for N in range(1,int(sys.argv[1])+1):
    I,T=build(N); far=set(I.far_pairs()); Ts=set(T)
    def Fpart(Tset):
        r=opt_cpsat(I,respect=Tset,free=far-Tset)
        return None if r is None else r[0]
    t0=time.time()
    full=Fpart(Ts); empty=Fpart(set()); drops=[Fpart(Ts-{e}) for e in T]
    print('N',N,'|T|',len(T),'partition F(T)-|H|=',full-len(I.H),'F(empty)-|H|=',empty-len(I.H),'F(T-e)-|H| values',sorted(set(d-len(I.H) for d in drops)),round(time.time()-t0,1),flush=True)
