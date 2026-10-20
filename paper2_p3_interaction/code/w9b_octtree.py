"""Tree of octahedral K4 units glued at vertices. unit k has a glue vertex g_k (local id) toward its parent; parent glue vertex pv_k (local id in the parent)."""
import itertools, sys, random, time
from collections import Counter
import numpy as np
from w9b_lib import *
from w9b_sigma import XT, max_core
from w9b_octchain import unit_members, perms, vid

def make_unit(rng):
    roles=tuple(rng.choice(perms) for _ in range(4))
    H=unit_members(roles)
    soft=set()
    for (a,b,c) in H: soft|={key(a,b),key(b,c),key(a,c)}
    hard=[p for p in itertools.combinations(range(6),2) if p not in soft]
    st=[rng.randrange(2) for _ in hard]
    E=set()
    for (a,b,c) in H: E.add(key(a,b)); E.add(key(b,c))
    for p,s in zip(hard,st):
        if s: E.add(p)
    return H,E,roles,tuple(st)

def build_tree(units, edges, cross):
    """units: list of (H,E); edges: list of (child_idx, child_local_vertex, parent_idx, parent_local_vertex); cross: list of ((i,x),(j,y)) hard edges"""
    gid={}; nxt=0
    parent_of={}
    for (c,cv,p,pv) in edges: parent_of[c]=(cv,p,pv)
    maps=[None]*len(units)
    for k in range(len(units)):
        mp={}
        for v in range(6):
            if k in parent_of and parent_of[k][0]==v:
                _,p,pv=parent_of[k]
                mp[v]=maps[p][pv]
            else:
                mp[v]=nxt; nxt+=1
        maps[k]=mp
    H=[];E=set()
    for k,(Hk,Ek) in enumerate(units):
        for (a,b,c) in Hk: H.append((maps[k][a],maps[k][b],maps[k][c]))
        for (u,v) in Ek: E.add(key(maps[k][u],maps[k][v]))
    for (i,x),(j,y) in cross:
        u=maps[i][x]; v=maps[j][y]
        if u!=v: E.add(key(u,v))
    return Inst(nxt,E,H),maps

def enum_sigma(I, limit=200000, drop_far=False):
    n=I.n; m=len(I.H)
    farset=set(I.far_pairs()) if drop_far else set()
    # order members by BFS adjacency for good pruning
    order=[0]; seen={0}
    adjm=[set() for _ in range(m)]
    for i in range(m):
        for v in I.H[i]:
            for j in I.mem[v]:
                if j!=i: adjm[i].add(j)
    while len(order)<m:
        cand=[j for j in range(m) if j not in seen and any(j in adjm[i] for i in order)]
        if not cand: cand=[j for j in range(m) if j not in seen]
        # choose the one with most neighbours in seen
        cand.sort(key=lambda j:-len(adjm[j]&seen)); order.append(cand[0]); seen.add(cand[0])
    pos={i:k for k,i in enumerate(order)}
    by_last=[[] for _ in range(m)]
    for x,y,z in itertools.combinations(range(n),3):
        ps=[key(x,y),key(y,z),key(x,z)]
        if any(p in farset for p in ps): continue
        ow=[I.owner.get(p) for p in ps]
        own={o[0] for o in ow if o}
        if not own:
            if sum(I.adj[p[0]][p[1]] for p in ps)==2: return []
            continue
        by_last[max(pos[o] for o in own)].append((ps,ow))
    st=[0]*m; out=[]
    def ok(k):
        for ps,ow in by_last[k]:
            s=0
            for p,o in zip(ps,ow):
                s+= (XT[o[1]][st[o[0]]] if o else I.adj[p[0]][p[1]])
            if s==2: return False
        return True
    def rec(k):
        if len(out)>limit: return
        if k==m: out.append(tuple(st)); return
        i=order[k]
        for s in (0,1,2):
            st[i]=s
            if ok(k): rec(k+1)
    rec(0)
    return out

def random_tree(rng, shape, kcross):
    """shape: 'chain' or 'star' or 'bin' ; returns Inst"""
    if shape=='star2':   # root + 2 children
        nu=3; edges=[(1,None,0,None),(2,None,0,None)]
    elif shape=='star3': nu=4; edges=[(1,None,0,None),(2,None,0,None),(3,None,0,None)]
    elif shape=='chain3': nu=3; edges=[(1,None,0,None),(2,None,1,None)]
    elif shape=='chain4': nu=4; edges=[(1,None,0,None),(2,None,1,None),(3,None,2,None)]
    units=[]; 
    for _ in range(nu):
        H,E,roles,st=make_unit(rng); units.append((H,E))
    # choose glue vertices
    used=[set() for _ in range(nu)]
    ed=[]
    for (c,_,p,_) in edges:
        cv=rng.randrange(6)
        while cv in used[c]: cv=rng.randrange(6)
        pv=rng.randrange(6)
        while pv in used[p]: pv=rng.randrange(6)
        used[c].add(cv); used[p].add(pv); ed.append((c,cv,p,pv))
    cross=[]
    for _ in range(kcross):
        (c,cv,p,pv)=rng.choice(ed)
        cross.append(((c,rng.randrange(6)),(p,rng.randrange(6))))
    I,maps=build_tree(units,ed,cross)
    return I

if __name__=='__main__':
    seed=int(sys.argv[1]); trials=int(sys.argv[2]); shape=sys.argv[3]
    rng=random.Random(seed); dist=Counter(); best=None; t0=time.time()
    for t in range(trials):
        I=random_tree(rng,shape,rng.choice([0,0,1,2,3]))
        if not I.valid() or I.cmax()>4: continue
        S=enum_sigma(I)
        if len(S)<2: continue
        Sg=np.array(S,dtype=np.int8)
        k1=max_core(Sg,singleton_only=True); k2=max_core(Sg)
        dist[(k1,k2)]+=1
        if best is None or (k1,k2)>best[0]: best=((k1,k2),(I.n,sorted(I.E),I.H))
    print(shape,'trials',trials,'dist',sorted(dist.items()),'time',round(time.time()-t0))
    print('best',best[0] if best else None)
