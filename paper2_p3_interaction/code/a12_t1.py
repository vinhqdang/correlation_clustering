"""Claim (1): 5-state vs 7-state model, checked against exact OPT; all optimal partitions' member states."""
import sys, itertools, random
from a12_lib import *
rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 1)
NT=int(sys.argv[2]) if len(sys.argv)>2 else 60
sets7=[m for m in range(1,8)]
def rgs(n):
    a=[0]*n
    def rec(i,m):
        if i==n: yield a; return
        for v in range(m+1):
            a[i]=v; yield from rec(i+1,max(m,v+1))
    yield from rec(0,0)
def pcost(n,E,lab):
    c=0
    for a in range(n):
        for b in range(a+1,n):
            if (K(a,b) in E)!=(lab[a]==lab[b]): c+=1
    return c
stat={'inst':0,'opt_ok5':0,'opt_ok7':0,'bad':0,'allopt_states':set(),'nsol':0,'min7_uses_nonvalid':0}
for it in range(NT):
    n=rng.randint(5,8); E=rand_graph(n,rng.choice([.3,.4,.5]),rng)
    H=rand_packing(n,E,rng,maximal=rng.random()<.5,drop=.4)
    H=H[:4]
    if not H: continue
    # packing restricted to <=4 members (still a packing, possibly not maximal)
    OPT=opt_dp(n,E); k=OPT-len(H)
    soft=soft_set(H); hard=[K(a,b) for a in range(n) for b in range(a+1,n) if K(a,b) not in soft]
    # (a) all optimal partitions: soft edit sets
    for lab in rgs(n):
        if pcost(n,E,lab)==OPT:
            F={K(a,b) for a in range(n) for b in range(a+1,n) if (K(a,b) in E)!=(lab[a]==lab[b])}
            stat['nsol']+=1
            tot=0
            for t in H:
                P=prs(t); m=sum(1<<i for i in range(3) if P[i] in F)
                assert m>0, 'unedited member'
                stat['allopt_states'].add(m); tot+=bin(m).count('1')
            W=F-soft; assert tot+len(W)==OPT
    # (b) models min
    def model(allowed_masks):
        best=10**9
        for kk in range(0,k+1):
            for W in itertools.combinations(hard,kk):
                EW=set(E)^set(W)
                for xs in itertools.product(allowed_masks,repeat=len(H)):
                    E2=set(EW); cost=kk
                    for t,m in zip(H,xs):
                        P=prs(t)
                        for i in range(3):
                            if m>>i&1: E2^={P[i]}
                        cost+=bin(m).count('1')
                    if cost<best and is_cluster(n,E2): best=cost
        return best
    m5=model([1,2,4,3,7]); m7=model(sets7)
    stat['inst']+=1
    stat['opt_ok5']+= (m5==OPT); stat['opt_ok7']+=(m7==OPT)
    if m5!=OPT or m7!=OPT: stat['bad']+=1; print('BAD',n,sorted(E),H,OPT,m5,m7)
print(stat)
