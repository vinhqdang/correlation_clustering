import sys, itertools, random
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from r1_lib import *
rng=random.Random(5)
cnt={};fail={}
def chk(k,c,i=None):
    cnt[k]=cnt.get(k,0)+1
    if not c: fail.setdefault(k,[]).append(i)
tight=0; tried=0
for it in range(1500):
    n=rng.choice([6,7,8]); E=rand_graph(n,rng.choice([0.3,0.45,0.6]),rng)
    if not p3s(n,E): continue
    H,used=rand_packing(n,E,rng,maximal=True,keep=rng.choice([0.3,1.0]))
    nH=len(H); 
    if not H: continue
    Q={v:set() for v in range(n)}
    for P in used:
        a,b=tuple(P); Q[a].add(b); Q[b].add(a)
    D=max(len(Q[v]) for v in range(n))
    adj={v:{u for u in range(n) if frozenset((u,v)) in E} for v in range(n)}
    hadj={v:{u for u in adj[v] if frozenset((u,v)) not in used} for v in range(n)}
    ap=allparts_cost(n,E); m=min(c for c,_ in ap); ell=m-nH
    VH=set(v for t in H for v in t); R=[v for v in range(n) if v not in VH]
    # shrink tightness at ell=1,  |X|=2
    seen=set()
    for r in R:
        if r in seen: continue
        X={r}; st=[r]
        while st:
            x=st.pop()
            for y in adj[x]:
                if y in R and y not in X: X.add(y); st.append(y)
        seen|=X
        if len(X)==2:
            x=max(X); keep=[v for v in range(n) if v!=x]; mp={v:i for i,v in enumerate(keep)}
            E2={frozenset((mp[a],mp[b])) for e in E for a,b in [tuple(e)] if a!=x and b!=x}
            o2=opt(n-1,E2)
            for L in (1,):
                tried+=1
                if (m<=nH+L)!=(o2<=nH+L): tight+=1
    B=2*D+3*max(ell,0)
    for c,p in ap:
        if c!=m: continue
        lab={}
        for i,b in enumerate(p):
            for v in b: lab[v]=i
        S={frozenset((a,b)) for a,b in itertools.combinations(range(n),2) if (frozenset((a,b)) in E)!=(lab[a]==lab[b])}
        U={e for e in S if e not in used}; u=len(U)
        exc=set(lab[v] for e in U for v in e)
        chk('alpha_count',len(exc)<=2*u)
        for i in exc:
            C=p[i]; chk('alpha_size',len(C)<=1+2*D+u,(n,E,H,p)); 
            chk('alpha_LB',all(len(hadj[v])<=B for v in C),(n,E,H,p))
            # connected
            Cs=set(C); comp={C[0]}; st=[C[0]]
            while st:
                x=st.pop()
                for y in adj[x]&Cs:
                    if y not in comp: comp.add(y); st.append(y)
            chk('conn',comp==Cs)
        for i,C in enumerate(p):
            if i in exc: continue
            Cs=set(C)
            chk('beta_noleave',all(hadj[v]<=Cs for v in C))
            # hard components inside
            comps=set()
            for v in C:
                comp={v}; st=[v]
                while st:
                    x=st.pop()
                    for y in hadj[x]:
                        if y not in comp: comp.add(y); st.append(y)
                comps.add(frozenset(comp))
            chk('beta_union',set().union(*comps)==Cs)
            if len(comps)>=2: chk('beta_2D',len(C)<=2*D,(n,E,H,p))
print({k:(cnt[k],'FAIL' if k in fail else 'ok') for k in cnt})
print('shrink |X|=2 at ell=1 tried',tried,'differs',tight)
