"""(a) star with m P3s through the hub: OPT, l, and I(H)=K_m;
(b) over random small graphs: in EVERY optimal partition each packed P3 receives >=1 edit,
    hence #unpacked edits <= OPT-|H| (the step behind the n^(2l) guess of unpacked edits)."""
import itertools,random
from check_csp_view import packing,pairs,partitions
def edits(lab,n,E):
    return {frozenset((a,b)) for a,b in itertools.combinations(range(n),2) if (lab[a]==lab[b])!=(frozenset((a,b)) in E)}
def opt_all(n,E):
    best=10**9;sols=[]
    for lab in partitions(n):
        S=edits(lab,n,E)
        if len(S)<best: best=len(S);sols=[S]
        elif len(S)==best: sols.append(S)
    return best,sols
for m in (2,3):
    n=2*m+1; E={frozenset((0,i)) for i in range(1,n)}
    H=[(2*i+1,0,2*i+2) for i in range(m)]
    o,_=opt_all(n,E); print("star leaves",2*m,"|H|",m,"OPT",o,"l",o-m)
rng=random.Random(5);cnt=0
for _ in range(300):
    n=rng.randint(4,7);p=rng.choice([.3,.5])
    E={frozenset(e) for e in itertools.combinations(range(n),2) if rng.random()<p}
    H=packing(n,E,rng)
    if not H: continue
    o,sols=opt_all(n,E)
    for S in sols:
        assert all(any(q in S for q in pairs(P)) for P in H)
        packed={q for P in H for q in pairs(P)}
        assert len(S-packed)<=o-len(H)
    cnt+=1
print("random instances checked",cnt,"ok")
