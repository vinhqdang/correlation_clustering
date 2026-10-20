"""Chains of octahedral K4 units glued at antipodal vertices (in-vertex of unit j+1 = out-vertex of unit j)."""
import itertools, pickle, sys, random, time
from collections import Counter
import numpy as np
from w9b_lib import *
from w9b_sigma import XT, max_core
pairsK=list(itertools.combinations(range(4),2)); vid={p:k for k,p in enumerate(pairsK)}
perms=list(itertools.permutations(range(3)))
ANTI=[(0,5),(1,4),(2,3)]    # antipodal pairs in the vertex numbering v_{01},v_{02},v_{03},v_{12},v_{13},v_{23} -> ids 0..5 : (01|23),(02|13),(03|12)
def unit_members(roles):
    H=[]
    for i in range(4):
        vs=[vid[tuple(sorted((i,j)))] for j in range(4) if j!=i]
        a,b,c=[vs[k] for k in roles[i]]
        H.append((a,b,c))
    return H
def triple_ok(adj_get, states_get, owner, tr):
    pass

class Unit:
    def __init__(self,roles,st,inv,outv):
        self.H=unit_members(roles)
        soft=set()
        for (a,b,c) in self.H: soft|={key(a,b),key(b,c),key(a,c)}
        self.hard=[p for p in itertools.combinations(range(6),2) if p not in soft]
        self.E=set()
        for (a,b,c) in self.H: self.E.add(key(a,b)); self.E.add(key(b,c))
        for p,s in zip(self.hard,st):
            if s: self.E.add(p)
        self.inv=inv; self.outv=outv
        I=Inst(6,self.E,self.H)
        self.valid=I.valid()
        self.S=None
        if self.valid:
            from w9b_sigma import sigma
            self.S=[tuple(r) for r in sigma(I,drop_far=False).tolist()]

def build_chain_instance(U,N,cross_eps):
    """vertices: unit j uses local ids 0..5 ; out-vertex of unit j identified with in-vertex of unit j+1.
    cross_eps: set of (j, x_local, y_local) hard edge pairs between unit j and unit j+1 (x in unit j, y in unit j+1) """
    gid={}; nxt=0
    maps=[]
    for j in range(N):
        mp={}
        for v in range(6):
            if j>0 and v==U.inv: mp[v]=maps[j-1][U.outv]
            else: mp[v]=nxt; nxt+=1
        maps.append(mp)
    H=[]; E=set()
    for j in range(N):
        for (a,b,c) in U.H: H.append((maps[j][a],maps[j][b],maps[j][c]))
        for (u,v) in U.E: E.add(key(maps[j][u],maps[j][v]))
    for (j,x,y) in cross_eps: E.add(key(maps[j][x],maps[j+1][y]))
    return Inst(nxt,E,H)

def chain_sigma(U,N,cross_eps):
    """DFS join using real triple tests between consecutive units"""
    I=build_chain_instance(U,N,cross_eps)
    if not I.valid(): return None,I
    # compat between consecutive units via sigma restricted: evaluate all triples touching both units
    from w9b_sigma import sigma
    m=len(U.H)
    rows=[]
    # simple: full DFS with partial constraint checks (triples whose owners all assigned)
    n=I.n
    triples=[]
    for x,y,z in itertools.combinations(range(n),3):
        ps=[key(x,y),key(y,z),key(x,z)]
        ow=[I.owner.get(p) for p in ps]
        own={o[0] for o in ow if o}
        if not own:
            if sum(I.adj[p[0]][p[1]] for p in ps)==2: return [],I
            continue
        triples.append((ps,ow,max(own)))
    by_last=[[] for _ in I.H]
    for ps,ow,l in triples: by_last[l].append((ps,ow))
    st=[0]*len(I.H); out=[]
    def ok(k):
        for ps,ow in by_last[k]:
            s=0
            for p,o in zip(ps,ow):
                s+= (XT[o[1]][st[o[0]]] if o else I.adj[p[0]][p[1]])
            if s==2: return False
        return True
    def rec(k):
        if k==len(I.H): out.append(tuple(st)); return
        for s in (0,1,2):
            st[k]=s
            if ok(k): rec(k+1)
    rec(0)
    return out,I

if __name__=='__main__':
    seed=int(sys.argv[1]); trials=int(sys.argv[2]); N=int(sys.argv[3])
    rng=random.Random(seed)
    best=Counter(); bestinfo={}
    t0=time.time()
    for t in range(trials):
        roles=tuple(rng.choice(perms) for _ in range(4))
        st=tuple(rng.randrange(2) for _ in range(3))
        inv,outv=rng.choice([(a,b) for a,b in itertools.permutations(range(6),2)])
        U=Unit(roles,st,inv,outv)
        if not U.valid or not U.S or len(U.S)<2: continue
        # members through in/out vertices must be 2 each; any vertex is in exactly 2 members
        ce=set()
        for _ in range(rng.choice([0,1,2,3])):
            ce.add((0,rng.randrange(6),rng.randrange(6))) if False else None
        # cross edges: between unit j and j+1, x in unit j, y in unit j+1 (not both the glue vertex)
        k=rng.choice([0,0,1,2,3,4])
        for _ in range(k):
            j=rng.randrange(N-1); x=rng.randrange(6); y=rng.randrange(6)
            if x==U.outv or y==U.inv: continue
            ce.add((j,x,y))
        S,I=chain_sigma(U,N,ce)
        if not S or len(S)<2: continue
        Sg=np.array(S,dtype=np.int8)
        k1=max_core(Sg,singleton_only=True); k2=max_core(Sg)
        key_=(k1,k2); best[key_]+=1
        if key_ not in bestinfo: bestinfo[key_]=(roles,st,inv,outv,sorted(ce),len(S))
    print('N',N,'trials',trials,'dist',sorted(best.items()),'time',round(time.time()-t0))
    top=max(best) if best else None
    print('top',top,bestinfo.get(top))
