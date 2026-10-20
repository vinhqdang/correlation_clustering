import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6c_nflib import *
import random, time
rng=random.Random(int(sys.argv[2])); n=int(sys.argv[1]); T=time.time(); best={}; tight=0; tot=0; hist={}
while time.time()-T<float(sys.argv[3]):
    m=rng.randint(n-1,int(2.3*n))
    g=nx.gnm_random_graph(n,m,seed=rng.randrange(10**9))
    if not nx.is_connected(g): continue
    A=np.array(nx.to_numpy_array(g),dtype=np.int8); tot+=1
    opt,lab,same,cost,P=brute_opt(A)
    optidx=np.nonzero(cost==opt)[0]
    Hs=packings_of_size(A,opt)
    if not Hs: continue
    tight+=1; rng.shuffle(Hs)
    for H in Hs[:20]:
        used=set()
        for t in H: used|=pairs_of(t)
        hg=nx.Graph(); hg.add_nodes_from(range(n))
        for i,j in P:
            if A[i][j] and (i,j) not in used: hg.add_edge(i,j)
        comps=[frozenset(c) for c in nx.connected_components(hg)]
        # closed components in each optimal partition
        fam=[]
        for i in optidx:
            l=lab[i]
            cl={}
            for v in range(n): cl.setdefault(l[v],set()).add(v)
            cls=[frozenset(b) for b in cl.values()]
            fam.append(frozenset(k for k,c in enumerate(comps) if c in cls))
        fam=list(set(fam))
        idxs=list(range(len(comps)))
        def isface(S): return any(S<=f for f in fam)
        nfs=[]
        for r in range(1,len(idxs)+1):
            for S in itertools.combinations(idxs,r):
                S=frozenset(S)
                if any(T_<=S for T_ in nfs): continue
                if not isface(S): nfs.append(S)
        mm=max((len(s) for s in nfs),default=0)
        hist[mm]=hist.get(mm,0)+1
        if mm>=3 and (n,mm) not in best:
            best[(n,mm)]=(sorted(g.edges()),H,[sorted(c) for c in comps],[sorted(s) for s in nfs if len(s)==mm]); print('FOUND',best[(n,mm)],flush=True)
print('n',n,'sampled',tot,'tight',tight,'maxnonface hist',hist,'found',list(best))
