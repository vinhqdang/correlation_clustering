import itertools, random
import numpy as np
_PART={}
def parts(n):
    if n in _PART: return _PART[n]
    res=[]
    def rec(i,lab,k):
        if i==n: res.append(tuple(lab)); return
        for b in range(k+1):
            lab.append(b); rec(i+1,lab,max(k,b+1) if b==k else k); lab.pop()
    rec(0,[],0)
    L=np.array(res,dtype=np.int8)
    prs=list(itertools.combinations(range(n),2))
    T=np.stack([(L[:,a]==L[:,b]) for a,b in prs],axis=1)
    _PART[n]=(L,prs,T); return _PART[n]
def rand_graph(n,p,rng):
    return set(frozenset(e) for e in itertools.combinations(range(n),2) if rng.random()<p)
def p3s(n,E):
    out=[]
    for y in range(n):
        for x,z in itertools.combinations([v for v in range(n) if v!=y],2):
            if frozenset((x,y)) in E and frozenset((y,z)) in E and frozenset((x,z)) not in E:
                out.append((x,y,z))
    return out
def pairs3(t):
    x,y,z=t
    return {frozenset((x,y)),frozenset((y,z)),frozenset((x,z))}
def rand_packing(n,E,rng,maximal=True,keep=1.0):
    ts=p3s(n,E); rng.shuffle(ts); used=set(); H=[]
    for t in ts:
        if rng.random()>keep: continue
        P=pairs3(t)
        if not (P&used): H.append(t); used|=P
    if maximal:
        for t in ts:
            P=pairs3(t)
            if not (P&used): H.append(t); used|=P
    return H,used
def costs(n,E):
    L,prs,T=parts(n)
    Ev=np.array([frozenset(p) in E for p in prs])
    return (T!=Ev).sum(1)
def blocks(lab):
    d={}
    for v,b in enumerate(lab): d.setdefault(int(b),[]).append(v)
    return list(d.values())
