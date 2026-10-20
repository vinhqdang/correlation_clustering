"""verify: every vertex triple with at least one non-far pair has its non-far-pair members inside a ball of radius <=3 of I(H)
(so all variables of the triple constraint live in one bag N^3[beta_t])"""
import sys, random, itertools, time
import networkx as nx
from collections import Counter
from w9b_lib import *

def rand_tree_instance(rng,m,cmax,extra):
    """forest/graph of members: member i shares a vertex with an earlier member with prob 0.8; random extra shared vertices; random hard edges"""
    H=[]; cnt={}; nv=0
    for i in range(m):
        vs=[]
        if i>0 and rng.random()<0.85:
            for _ in range(20):
                j=rng.randrange(max(0,i-4),i)
                v=rng.choice(H[j])
                if cnt.get(v,0)<cmax: vs.append(v); break
        while len(vs)<3:
            if len(vs)>=1 and rng.random()<extra and i>0:
                j=rng.randrange(i); v=rng.choice(H[j])
                if v not in vs and cnt.get(v,0)<cmax: vs.append(v); continue
            vs.append(nv); nv+=1
        t=tuple(rng.sample(vs,3)); H.append(t)
        for v in t: cnt[v]=cnt.get(v,0)+1
    # reject duplicate pairs
    seen=set()
    for (a,b,c) in H:
        for p in (key(a,b),key(b,c),key(a,c)):
            if p in seen: return None
            seen.add(p)
    n=nv+rng.randint(0,3)
    E=set()
    for (a,b,c) in H: E.add(key(a,b)); E.add(key(b,c))
    p=rng.choice([0.05,0.15,0.3])
    for u,v in itertools.combinations(range(n),2):
        if key(u,v) in seen: continue
        if rng.random()<p: E.add((u,v))
    return Inst(n,E,H)

def check(I):
    g=I.IH(); d=dict(nx.all_pairs_shortest_path_length(g))
    far=set(I.far_pairs())
    worst=0
    m=len(I.H)
    for x,y,z in itertools.combinations(range(I.n),3):
        ps=[key(x,y),key(y,z),key(x,z)]
        nonfar=[p for p in ps if p not in far]
        if not nonfar: continue
        verts={v for p in nonfar for v in p}
        W=set()
        for v in verts: W|=set(I.mem[v])
        if not W: continue
        r=min(max(d[c].get(w,99) for w in W) for c in range(m) if c in d and (not W or True))
        worst=max(worst,r)
    return worst

if __name__=='__main__':
    seed=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); cmax=int(sys.argv[4])
    rng=random.Random(seed); dist=Counter(); t0=time.time()
    for _ in range(N):
        I=rand_tree_instance(rng,m,cmax,rng.choice([0.0,0.2,0.4]))
        if I is None or not I.valid(): continue
        dist[check(I)]+=1
    print('m',m,'c<=',cmax,'max over triples of the min-centre radius, distribution:',sorted(dist.items()),round(time.time()-t0),'s')
