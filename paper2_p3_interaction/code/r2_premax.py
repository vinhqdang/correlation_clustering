# audit9: brute-force of Lemma 'Reduction preserves pair-maximality' (lem:premax) and counting Lemma 5.6 (lem:count)
import sys, random, itertools
from r2_lib import *
seed=int(sys.argv[1]); N=int(sys.argv[2]); rng=random.Random(seed)
def comps_of(H):
    par=list(range(len(H)))
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i,j in itertools.combinations(range(len(H)),2):
        if set(H[i])&set(H[j]): par[f(i)]=f(j)
    d={}
    for i in range(len(H)): d.setdefault(f(i),[]).append(i)
    return list(d.values())
def sub_opts(n,E,W):
    W=sorted(W); m=len(W); idx={v:i for i,v in enumerate(W)}
    E2=set(frozenset((idx[a],idx[b])) for e in E for a,b in [tuple(e)] if a in idx and b in idx)
    cs=costs(m,E2); opt=int(cs.min()); L,prs,T=parts(m)
    sols=[ [ [W[i] for i in b] for b in blocks(L[k])] for k in np.nonzero(cs==opt)[0]]
    return opt,sols
def boundary_ok(n,E,W,bl):
    Wset=set(W)
    for y in range(n):
        if y in Wset: continue
        A={x for x in W if frozenset((x,y)) in E}
        if A and not any(A==set(b) for b in bl): return False
    return True
def apply(E,bl,W):
    E2=set(E)
    for a,b in itertools.combinations(sorted(W),2):
        tog=any(a in c and b in c for c in bl)
        pr=frozenset((a,b))
        if (pr in E)!=tog: E2^={pr}
    return E2
st=dict(red=0,irr_tight=0,premax_checked=0,count_checked=0,viol=0,opt_ok=0)
for it in range(N):
    n=rng.randint(6,9); p=rng.choice([0.2,0.3,0.4,0.5])
    E=rand_graph(n,p,rng)
    pm=rng.random()<0.6
    H,used=rand_packing(n,E,rng,maximal=pm,keep=rng.choice([1.0,0.7]))
    if len(H)<1: continue
    cl=comps_of(H); info=[]
    for c in cl:
        W=sorted({v for i in c for v in H[i]})
        tau,sols=sub_opts(n,E,W)
        good=[s for s in sols if boundary_ok(n,E,W,s)]
        info.append((c,W,tau,good))
    cs=costs(n,E); opt=int(cs.min())
    # lemma count: all partitions
    sumtau=sum(t for _,_,t,_ in info); J=sum(1 for _,_,_,g in info if not g)
    L,prs,T=parts(n)
    if (cs < sumtau+J/2-1e-9).any():
        st['viol']+=1; print("COUNT VIOL",seed,n,sorted(map(sorted,E)),H)
    st['count_checked']+=1
    if opt<sumtau: st['viol']+=1; print("OPT<sumtau")
    for c,W,tau,good in info:
        if tau==len(c) and not good: st['irr_tight']+=1
        if tau==len(c) and good and pm:
            # premax lemma, with any good S_F
            for S in good:
                st['red']+=1
                E2=apply(E,S,W); H2=[H[i] for i in range(len(H)) if i not in c]
                ok=True
                used2=set()
                for t in H2:
                    x,y,z=t
                    if not (frozenset((x,y)) in E2 and frozenset((y,z)) in E2 and frozenset((x,z)) not in E2): ok=False
                    P=pairs3(t)
                    if P&used2: ok=False
                    used2|=P
                if not ok: st['viol']+=1; print("not packing",seed)
                for t in p3s(n,E2):
                    if not (pairs3(t)&used2): st['viol']+=1; print("PREMAX VIOL",seed,n,sorted(map(sorted,E)),H,c,S); break
                st['premax_checked']+=1
                o2=int(costs(n,E2).min())
                if o2+len(c)==opt: st['opt_ok']+=1
                else: st['viol']+=1; print("OPT eq viol",seed)
print(seed,st)
