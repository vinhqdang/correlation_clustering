import random, sys, itertools
from collections import Counter
from w8b_lib import Inst
from w8b_ilp import ILP
from w8b_branch import gen_chain
def exact(o):
    I=o.I
    cnt=o.counted; o.counted=[True]*len(cnt); 
    v=o.F({}, []); o.counted=cnt; return v
seed=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); rng=random.Random(seed)
mc=Counter(); vc=Counter(); bad=0; tot=0; maxv=0
for _ in range(N):
    I=gen_chain(m,rng.randint(2,2+m),rng)
    if I is None: continue
    o=ILP(I); opt=exact(o); ex=opt-len(I.H)
    tot+=1
    for l in range(0,min(ex+1,3)):
        r,st=o.run(l)
        if r!=(ex<=l): bad+=1
        for c in st['cores']: mc[c]+=1
print("m",m,"inst",tot,"bad",bad,"core size dist",sorted(mc.items()))
