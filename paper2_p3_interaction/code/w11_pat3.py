import sys, time, itertools
from collections import Counter
from w11_lib import *
from w11_pat import build, NAMES
N=int(sys.argv[1]); l=int(sys.argv[2])
allsub=[''.join(c) for k in range(0,7) for c in itertools.combinations(NAMES,k)]
def ok(kind,S):
    I=build(3,(S,) if kind=='mid' else (),(S,) if kind=='start' else (),(S,) if kind=='end' else ())
    return I.valid() and I.pair_maximal()
val={k:[S for S in allsub if S and ok(k,S)] for k in ('mid','start','end')}
res=Counter(); best=[]; t0=time.time()
# two start guards, two end guards, one mid pendant pattern (more guards -> larger excess)
for S in val['mid']+['']:
    for st in itertools.combinations_with_replacement(val['start'],2):
        for en in itertools.combinations_with_replacement(val['end'],2):
            mid=(S,) if S else ()
            I=build(N,mid,st,en)
            if not (I.valid() and I.pair_maximal()): continue
            sz,T=min_certificate(I,l,max_size=30,workers=2)
            res[sz]+=1
            if sz is not None: best.append((sz,mid,st,en))
            if time.time()-t0>1500: break
best.sort(reverse=True)
print('N',N,'l',l,res,round(time.time()-t0))
for b in best[:6]: print(b)
