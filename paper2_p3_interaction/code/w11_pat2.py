import sys, time, itertools, random
from collections import Counter
from w11_lib import *
from w11_pat import build, NAMES
N = int(sys.argv[1])
allsub = [''.join(c) for k in range(0, 7) for c in itertools.combinations(NAMES, k)]
# valid single-pendant patterns (pm) tested on a middle-ish position in a 3-chain with the pendant on block 1 (mid) / block 0 (start) / block N-1 (end)
def ok(kind, S):
    I = build(3, (S,) if kind=='mid' else (), (S,) if kind=='start' else (), (S,) if kind=='end' else ())
    return I.valid() and I.pair_maximal()
val = {k: [S for S in allsub if S and ok(k, S)] for k in ('mid','start','end')}
print({k: len(v) for k,v in val.items()}); print(val['mid']); print(val['start']); print(val['end'])
t0=time.time()
res = Counter(); best=[]
for S in val['mid'] + ['']:
    for st in val['start']:
        for en in val['end']:
            mid = (S,) if S else ()
            I = build(N, mid, (st,), (en,))
            if not (I.valid() and I.pair_maximal()): continue
            sz, T = min_certificate(I, 0, max_size=16)
            res[sz] += 1
            if sz is not None: best.append((sz, mid, st, en))
best.sort(reverse=True)
print('N',N,res, round(time.time()-t0),'s')
for b in best[:20]: print(b)
