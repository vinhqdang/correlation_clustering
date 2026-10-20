import random, itertools, sys
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2])
def run():
    m=rng.randint(4,10); par=[-1]+[rng.randrange(v) if rng.random()<.9 else -1 for v in range(1,m)]
    depth=[0]*m
    for v in range(m): depth[v]=0 if par[v]<0 else depth[par[v]]+1
    def anc(v):
        a=[v]
        while par[a[-1]]>=0: a.append(par[a[-1]])
        return a
    def path(u,v):
        au,av=anc(u),anc(v)
        if au[-1]!=av[-1]: return None
        sv=set(av); l=next(x for x in au if x in sv)
        return set(au[:au.index(l)+1])|set(av[:av.index(l)+1]),l
    confl=[]
    for _ in range(rng.randint(1,7)):
        u,v=rng.randrange(m),rng.randrange(m); p=path(u,v)
        if p: confl.append((u,v,p[0],p[1]))
    d=rng.choice([1,2,3]); deg=[0]*m; links=[]
    for _ in range(rng.randint(0,m)):
        u,v=rng.randrange(m),rng.randrange(m)
        if u!=v and deg[u]<d and deg[v]<d and (u,v) not in links and (v,u) not in links: links.append((u,v)); deg[u]+=1; deg[v]+=1
    ops=[('b',v) for v in range(m)]+[('t',l) for l in links]
    def hits(op,c):
        u,v,P,l=c
        return (op[1] in P) if op[0]=='b' else (u in op[1] or v in op[1])
    def feasible(S): return all(any(hits(o,c) for o in S) for c in confl)
    opt=None
    for k in range(0,len(confl)+1):
        if any(feasible(S) for S in itertools.combinations(ops,k)): opt=k;break
    # search tree
    calls=[0]
    def st(S,L):
        calls[0]+=1
        un=[c for c in confl if not any(hits(o,c) for o in S)]
        if not un: return True
        if L==0: return False
        Q=max(un,key=lambda c:depth[c[3]])
        if st(S+[('b',Q[3])],L-1): return True
        for e in (Q[0],Q[1]):
            for l in links:
                if e in l and st(S+[('t',l)],L-1): return True
        return False
    for L in range(0,len(confl)+1):
        if st([],L): break
    return opt,L,d
bad=0; tun=0
for _ in range(N):
    o,L,d=run()
    if o!=L: bad+=1; print('MISMATCH',o,L)
print('RM2 searchtree instances',N,'mismatch',bad)
