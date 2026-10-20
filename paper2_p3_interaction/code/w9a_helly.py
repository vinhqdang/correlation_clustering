"""Helly-3 test: forest CSP on ZE states {L,R,M}; edge relations drawn from the 144 configurations of Lemma 1 (role pair x 2x2 cross bits);
arbitrary unary lists on arbitrary members.  Every minimal infeasible sub-family of lists must have <= 3 lists."""
import random, itertools, sys
ST={'L':(0,1,0),'R':(1,0,0),'M':(1,1,1)}
def view(s,role):
    ab,bc,ac=ST[s]
    return {'a':(ab,ac),'b':(ab,bc),'c':(bc,ac)}[role]
def compat(vu,vv,eps):
    for i in range(2):
        for j in range(2):
            if eps[i][j]:
                if vu[i]!=vv[j]: return False
            else:
                if vu[i] and vv[j]: return False
    return True
STS=['L','R','M']
RELS=set()
for ra in 'abc':
    for rb in 'abc':
        for bits in itertools.product((0,1),repeat=4):
            eps=((bits[0],bits[1]),(bits[2],bits[3]))
            R=frozenset((s,t) for s in STS for t in STS if compat(view(s,ra),view(t,rb),eps))
            RELS.add(R)
RELS=sorted(RELS,key=lambda r:(len(r),sorted(r)))
print('distinct ZE relations',len(RELS),'sizes',sorted(set(len(r) for r in RELS)))
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); m=int(sys.argv[3])
def feasible(par,rel,lists):
    # tree CSP, nodes 0..m-1, par[v]<v or -1; domain STS
    dom=[set(l) for l in lists]
    # DP bottom-up (nodes in reverse index order)
    ok=[set(d) for d in dom]
    for v in range(m-1,0,-1):
        p=par[v]
        if p<0: continue
        keep=set()
        for s in ok[p]:
            if any((s,t) in rel[v] for t in ok[v]): keep.add(s)
        ok[p]=keep
    return all(ok[v] for v in range(m) if par[v]<0)
mx=0; hist={}
for it in range(N):
    par=[-1]+[rng.randrange(v) if rng.random()<0.9 else -1 for v in range(1,m)]
    rel=[None]+[rng.choice(RELS) for _ in range(1,m)]
    # random list family: each member gets 0-2 lists, each a nonempty proper subset
    fam=[]
    for v in range(m):
        for _ in range(rng.choice([0,1,1,2])):
            fam.append((v,frozenset(rng.choice([('L',),('R',),('M',),('L','R'),('L','M'),('R','M')]))))
    full=[set(STS)]*m
    def feas(F):
        lists=[set(STS) for _ in range(m)]
        for (v,l) in F: lists[v]&=set(l)
        return feasible(par,rel,lists)
    if feas(fam): continue
    # minimal infeasible subfamily by deletion with random orders
    for _ in range(3):
        F=list(range(len(fam))); rng.shuffle(F)
        for x in list(F):
            G=[y for y in F if y!=x]
            if not feas([fam[i] for i in G]): F=G
        mx=max(mx,len(F)); hist[len(F)]=hist.get(len(F),0)+1
print('max minimal infeasible list family size',mx,hist)
