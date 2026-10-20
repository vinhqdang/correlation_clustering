import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6c_nflib import *
import random, time
rng=random.Random(int(sys.argv[2]) if len(sys.argv)>2 else 1)
n=int(sys.argv[1]); T=time.time(); best={}; tight=0; tot=0
while time.time()-T<float(sys.argv[3]) :
    m=rng.randint(n-1,int(1.9*n)); 
    g=nx.gnm_random_graph(n,m,seed=rng.randrange(10**9))
    if not nx.is_connected(g): continue
    A=np.array(nx.to_numpy_array(g),dtype=np.int8); tot+=1
    opt,lab,same,cost,P=brute_opt(A)
    optidx=np.nonzero(cost==opt)[0]
    sing=[]
    for i in optidx:
        l=lab[i]; sizes={}
        for v in range(n): sizes.setdefault(l[v],[]).append(v)
        sing.append(frozenset(v for b in sizes.values() if len(b)==1 for v in b))
    sing=list(set(sing))
    Hs=packings_of_size(A,opt)
    if not Hs: continue
    tight+=1
    rng.shuffle(Hs)
    for H in Hs[:30]:
        used=set()
        for t in H: used|=pairs_of(t)
        anch=[v for v in range(n) if all(key(v,u) in used for u in range(n) if A[v][u])]
        if len(anch)<3: continue
        nf=minimal_nonfaces(sing,anch)
        mm=max((len(s) for s in nf),default=0)
        if mm>=3 and (n,mm) not in best:
            best[(n,mm)]=(sorted(g.edges()),H,anch,[sorted(s) for s in nf if len(s)==mm]); print('FOUND',best[(n,mm)],flush=True)
print('n',n,'sampled',tot,'tight',tight,'found',list(best))
