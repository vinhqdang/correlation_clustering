import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import evaluate, e_relaxed
seed=int(sys.argv[1]); N=int(sys.argv[2]); mmax=int(sys.argv[3])
rng=random.Random(seed)
stats={}
cnt=0; viol=0
while cnt<N:
    m=rng.randint(3,mmax)
    inst=gen_forest(rng,m,extra=rng.randint(0,3),p=rng.choice([0.1,0.2,0.3,0.4]),share=rng.choice([0.8,1.0]))
    if not nx.is_forest(inst.Ig()) or inst.soft_triangles(): continue
    rc,res=evaluate(inst.n,sorted(inst.E),inst.H)
    if rc!=0 or res[0]>=100000 or res[0]==0: continue
    full=(1<<m)-1
    g=res[0]
    cnt+=1
    red=lambda Z: g-min(e_relaxed(full,Z),g+100)   # may be inf -> treat
    single={}
    for A in range(m):
        v=e_relaxed(full,1<<A)
        single[A]=g-v if v<100000 else 0
    for k in (2,3):
        for Z in itertools.combinations(range(m),k):
            zm=sum(1<<a for a in Z)
            v=e_relaxed(full,zm)
            r=g-v
            s=sum(single[a] for a in Z)
            key_=(k,r-s)
            stats[key_]=stats.get(key_,0)+1
            if r>s and viol<10:
                viol+=1; print('SUPERADD',Z,'red',r,'sum singles',s,'g',g,inst.n,sorted(inst.E),inst.H)
print(sorted(stats.items()))
