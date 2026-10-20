# audit9: brute-force check of Lemma half / EX / NF / CAND (inserted pairs within G-distance 2, inside L_B) on ALL optimal partitions
import sys, random, itertools
from r2_lib import *
seed=int(sys.argv[1]); N=int(sys.argv[2]); rng=random.Random(seed)
GAM=[0]
viol=0; tested=0; ninsert=0; ndel=0; maxd=0; nU=0
for it in range(N):
    n=rng.randint(5,9); p=rng.choice([0.25,0.35,0.5,0.65,0.8,0.9])
    E=rand_graph(n,p,rng)
    H,used=rand_packing(n,E,rng,maximal=True)
    if not H: continue
    cs=costs(n,E); opt=int(cs.min()); ell=opt-len(H)
    Q={v:set() for v in range(n)}
    for t in H:
        for pr in pairs3(t):
            a,b=tuple(pr); Q[a].add(b); Q[b].add(a)
    D=max(len(Q[v]) for v in range(n))
    soft=used
    Nh={v:{u for u in range(n) if u!=v and frozenset((u,v)) in E and frozenset((u,v)) not in soft} for v in range(n)}
    B=2*D+3*ell; Z=B+D+1
    LB={v for v in range(n) if len(Nh[v])<=B}
    L,prs,T=parts(n)
    for idx in np.nonzero(cs==opt)[0]:
        lab=L[idx]; bl=blocks(lab); tested+=1
        cl={v:int(lab[v]) for v in range(n)}
        U=[(a,b) for (a,b) in prs if frozenset((a,b)) not in soft and ((frozenset((a,b)) in E)!=(cl[a]==cl[b]))]
        u=len(U); nU+=u
        def bad(msg):
            global viol; viol+=1; print("VIOL",msg,seed,n,sorted(map(sorted,E)),H,bl)
        if u>ell: bad("u>ell")
        dele={v:0 for v in range(n)}; ins={v:0 for v in range(n)}
        for a,b in U:
            tgt=ins if frozenset((a,b)) not in E else dele
            tgt[a]+=1; tgt[b]+=1
        exc=set()
        for C in bl:
            Cs=set(C)
            for v in C:
                dC=len([x for x in C if x!=v and frozenset((v,x)) in E])
                if 2*dC<len(C)-1: bad("half a")
            for a,b in itertools.combinations(C,2):
                if frozenset((a,b)) not in E and not any(z!=a and z!=b and frozenset((a,z)) in E and frozenset((z,b)) in E for z in C): bad("half b")
        for a,b in U:
            C=set(bl[0]) if False else None
            C=[c for c in bl if a in c][0]
            if frozenset((a,b)) in E: # deleted
                ndel+=1
                for x,y in ((a,b),(b,a)):
                    if not (len(Nh[x])<=D+dele[x]+dele[y]-1): bad("EX a del")
                    if len(Nh[x])>D+u: bad("EX a del<=D+u")
            else:
                ninsert+=1
                for x,y in ((a,b),(b,a)):
                    if not (len(Nh[x])<=D+ins[y]-1+dele[x]): bad("EX a ins")
                    if len(Nh[x])>D+u: bad("EX a ins<=D+u")
                # distance 2 in G inside L_B
                if a not in LB or b not in LB: bad("endpoint not in LB")
                ok=any(frozenset((a,z)) in E and frozenset((z,b)) in E and z in LB for z in C)
                if not ok: bad("no witness in LB")
                maxd=max(maxd,2)
            for x in (a,b):
                Cx=[c for c in bl if x in c][0]
                exc.add(tuple(Cx))
                if len(Cx)>1+D+len(Nh[x])+ins[x]: bad("EX b")
                if len(Cx)>1+2*D+u: bad("EX b2")
                for z in Cx:
                    if len(Nh[z])>len(Cx)-1+dele[z] or len(Nh[z])>2*D+2*u: bad("EX c")
                    if z not in LB: bad("exc cluster not in LB")
        if len(exc)>2*u: bad("#exc")
        # hard components
        comp={}
        hv=[v for v in range(n)]
        par=list(range(n))
        def f(x):
            while par[x]!=x: par[x]=par[par[x]]; x=par[x]
            return x
        for a,b in itertools.combinations(range(n),2):
            if frozenset((a,b)) in E and frozenset((a,b)) not in soft: par[f(a)]=f(b)
        for C in bl:
            if tuple(C) in exc: continue
            comps={f(v) for v in C}
            for a,b in itertools.combinations(C,2):
                pr=frozenset((a,b))
                if pr not in soft and ((pr in E)!=True): bad("beta: hard nonedge inside nonexc")
            for v in C:
                for x in Nh[v]:
                    if x not in C: bad("beta: hard edge leaves nonexc")
            if len(comps)>=2 and len(C)>2*D: bad("beta size")
        for K in set(f(v) for v in range(n) if len(Nh[v])>B):
            nU+=0; GAM[0]+=1
            Kv=[v for v in range(n) if f(v)==K]
            if sorted(Kv) not in [sorted(c) for c in bl]: bad("gamma")
print("seed",seed,"tested optima",tested,"insert",ninsert,"del",ndel,"viol",viol,"gamma checks",GAM[0])
