import itertools, pickle
from w9b_lib import *
from collections import Counter
# killer config: A=(0,1,2), B=(3,4,1), C=(0,3,5), hard (ab,ac,bc,(a,r),(b,p),(c,q)) = (0,0,None,0,1,0)
# vertices p=0,q=1,a=2,r=3,b=4,c=5 ; add D on {4,5,6}
def build(A,B,C,D,hardstat):
    """hardstat: dict pair->0/1 for all pairs not soft"""
    mem=[A,B,C,D]
    E=set()
    for (x,y,z) in mem:
        E.add(key(x,y)); E.add(key(y,z))
    soft=set()
    for (x,y,z) in mem:
        soft|={key(x,y),key(y,z),key(x,z)}
    n=7
    for p,s in hardstat.items():
        if s: E.add(key(*p))
    return Inst(n,E,mem)
cnt=Counter(); found=[]
hardpairs=[(2,4),(2,5),(4,5),(2,3),(4,0),(5,1)]  # as in scope enumeration (ab,ac,bc,(a,r),(b,p),(c,q))
A=(0,1,2);B=(3,4,1);C=(0,3,5)
for Dr in itertools.permutations((4,5,6)):
    D=Dr
    # soft pairs of D
    soft={key(*p) for p in itertools.combinations(D,2)}
    for A_ in [A]:
        base=[(2,4),(2,5),(2,3),(4,0),(5,1)]  # hard pairs from scope (excluding bc which is soft in D)
        extra=[(6,0),(6,1),(6,2),(6,3)]   # d with p,q,a,r
        allp=base+extra
        for st in itertools.product((0,1),repeat=len(allp)):
            hs=dict(zip(allp,st))
            # required statuses of the killer config: ab=0,ac=0,(a,r)=0,(b,p)=1,(c,q)=0 -> pairs (2,4)=0,(2,5)=0,(2,3)=0,(4,0)=1,(5,1)=0
            if hs[(2,4)]!=0 or hs[(2,5)]!=0 or hs[(2,3)]!=0 or hs[(4,0)]!=1 or hs[(5,1)]!=0: continue
            hh=dict(hs)
            # also pairs (6, other) all given; check all 21 pairs covered
            try:
                I=build(A,B,C,D,hh)
            except Exception as e: continue
            if not I.valid(): continue
            S=ze_states(I)
            cnt[len(S)]+=1
            found.append((D,st,S))
print(cnt)
pickle.dump(found,open('/tmp/w9b_diamond.pkl','wb'))
