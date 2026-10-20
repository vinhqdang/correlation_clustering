import sys, time, itertools, random
from collections import Counter
from w11_lib import *
from w11_pat import build, NAMES
l=int(sys.argv[1]); Ns=[int(x) for x in sys.argv[2].split(',')]; nsample=int(sys.argv[3])
allsub=[''.join(c) for k in range(0,7) for c in itertools.combinations(NAMES,k)]
def ok(kind,S):
    I=build(3,(S,) if kind=='mid' else (),(S,) if kind=='start' else (),(S,) if kind=='end' else ())
    return I.valid() and I.pair_maximal()
val={k:[S for S in allsub if S and ok(k,S)] for k in ('mid','start','end')}
rnd=random.Random(5)
combos=[]
for S in val['mid']+['']:
    for st in itertools.combinations_with_replacement(val['start'],2):
        for en in itertools.combinations_with_replacement(val['end'],2):
            combos.append(((S,) if S else (),st,en))
combos=[c for c in combos if (lambda I: I.valid() and I.pair_maximal())(build(3,*c))]
print('valid combos',len(combos))
rnd.shuffle(combos)
t0=time.time(); growth=Counter(); ex=[]
for mid,st,en in combos[:nsample]:
    sizes=[]
    for N in Ns:
        I=build(N,mid,st,en)
        if not (I.valid() and I.pair_maximal()): sizes=None; break
        sz,T=min_certificate(I,l,max_size=40,workers=2); sizes.append(sz)
        if sz is None: break
    if sizes is None: continue
    growth[tuple(sizes)]+=1
    if None not in sizes and sizes[-1]>sizes[0]: ex.append((sizes,mid,st,en))
print('l',l,'Ns',Ns,round(time.time()-t0)); print(sorted(growth.items(),key=str)); print('growth examples',ex[:5])
