import sys, itertools, random; sys.path.insert(0,'.')
from r1_lib import *
rng=random.Random(7); tested=0; bad=0; found=0
for it in range(1500):
    n=rng.choice([6,7,8]); E=rand_graph(n,rng.choice([0.2,0.3,0.4,0.5]),rng)
    if not p3s(n,E): continue
    H,used=rand_packing(n,E,rng,maximal=True,keep=rng.choice([0.4,1.0]))
    if not H: continue
    own={}
    for i,t in enumerate(H):
        for P in pairs3(t): own[P]=i
    Q={v:set() for v in range(n)}
    for P in used:
        a,b=tuple(P); Q[a].add(b); Q[b].add(a)
    D=max(len(Q[v]) for v in range(n))
    adj={v:{u for u in range(n) if frozenset((u,v)) in E} for v in range(n)}
    ap=allparts_cost(n,E); m=min(c for c,_ in ap)
    for c,p in ap:
        if c!=m: continue
        lab={}
        for i,b in enumerate(p):
            for v in b: lab[v]=i
        S={frozenset((a,b)) for a,b in itertools.combinations(range(n),2) if (frozenset((a,b)) in E)!=(lab[a]==lab[b])}
        U={e for e in S if e not in used}
        par=list(range(n))
        def f(x):
            while par[x]!=x: par[x]=par[par[x]]; x=par[x]
            return x
        for a in range(n):
            for b in adj[a]:
                e=frozenset((a,b))
                if a<b and e not in used and e not in U: par[f(a)]=f(b)
        for cl in p:
            pcs={}
            for v in cl: pcs.setdefault(f(v),[]).append(v)
            for A,B in itertools.combinations(list(pcs.values()),2):
                if any(frozenset((a,b)) in used for a in A for b in B): continue
                found+=1
                for a in A:
                    cnt=0
                    for b in B:
                        tested+=1
                        e=frozenset((a,b))
                        ok = (e in U and e not in E)
                        zs=[z for z in cl if z not in A and z not in B and frozenset((a,z)) in used and frozenset((z,b)) in used and own[frozenset((a,z))]!=own[frozenset((z,b))]]
                        ok = ok and bool(zs)
                        if not ok: bad+=1
                    # partners bound
                    cand={b for z in Q[a] for b in Q[z]} 
                print
print('softfree merges',found,'pair checks',tested,'violations',bad)
