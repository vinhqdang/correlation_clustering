"""Abstract tree-CSP calculus (forest I(H), c<=2): propagate (locked, unlocked) forced-state sets along chains of members,
with ALL constraints that an edge between adjacent members imposes (binary conj + the unary side constraints caused by the 4 cross pairs eps_{pr}),
then look for a hub (3 branches at the three vertices of mu) whose locks form a core of size 3."""
import itertools, io, contextlib, sys
with contextlib.redirect_stdout(io.StringIO()):
    from w7_rel import X, ok3, S
D=frozenset(S)
def pair(u,v): return ''.join(sorted((u,v)))
def edgecfg(qm,qn,eps):
    """parent member P (role qm at which child attaches), child C (role qn). eps indexed by (p,r), p in P\{qm}, r in C\{qn}.
    returns R subset of S x S (P state, C state), LamP, LamC."""
    om=[x for x in 'abc' if x!=qm]; on=[x for x in 'abc' if x!=qn]
    R={(s,t) for s in S for t in S}
    E={}
    k=0
    for p in om:
        for r in on:
            E[(p,r)]=eps[k]; k+=1
    for p in om:
        for r in on:
            R={(s,t) for (s,t) in R if ok3(X[pair(p,qm)][s],X[pair(qn,r)][t],E[(p,r)])}
    LamP=set(S)
    for r in on:   # triple (p1,p2,r), soft pair (p1,p2) of P
        LamP={s for s in LamP if ok3(X[pair(om[0],om[1])][s],E[(om[0],r)],E[(om[1],r)])}
    LamC=set(S)
    for p in om:   # triple (p, r1, r2), soft pair (r1,r2) of C
        LamC={t for t in LamC if ok3(X[pair(on[0],on[1])][t],E[(p,on[0])],E[(p,on[1])])}
    return R,frozenset(LamP),frozenset(LamC),E
cfgs={}
for qm in 'abc':
    for qn in 'abc':
        for eps in itertools.product((0,1),repeat=4):
            cfgs[(qm,qn,eps)]=edgecfg(qm,qn,eps)
def img(R,T): return frozenset(s for (s,t) in R if t in T)
# anchor leaf: role toward parent rho; anchor at end e in {a,c}, e!=rho; forced state: anchor at a -> L, at c -> R
# base[rho] = set of (TL,TU,anchorvertex) for a leaf member that has role rho toward parent
def leaf_states(rho):
    out=[]
    for e in 'ac':
        if e==rho: continue
        out.append((e, frozenset('L' if e=='a' else 'R')))
    return out
# branch = chain: we compute for a node with role rho toward parent the set of achievable (TL,TU) pairs (TU = unlocked).
from functools import lru_cache
@lru_cache(None)
def node(rho,depth):
    """achievable (TL,TU,desc) for a member with role rho toward its parent, subtree depth <= depth (leaf = depth 1)"""
    res={}
    # leaf
    for e,forced in leaf_states(rho):
        res.setdefault((forced,D),('leaf',rho,e))
    if depth>1:
        # internal chain node with one child attached at role sigma != rho
        for sigma in 'abc':
            if sigma==rho: continue
            for crho in 'abc':
                for (TLc,TUc),desc in node_items(crho,depth-1):
                    for eps in itertools.product((0,1),repeat=4):
                        R,LP,LC,E=cfgs[(sigma,crho,eps)]
                        # if child is a leaf, anchor vertex pairs must have eps 0
                        if desc[0]=='leaf':
                            e=desc[2]; on=[x for x in 'abc' if x!=crho]
                            if any(E[(p,e)]==1 for p in [x for x in 'abc' if x!=sigma]): continue
                        TL=LP & img(R,TLc&LC); TU=LP & img(R,TUc&LC)
                        if not TU: continue
                        key=(TL,TU)
                        res.setdefault(key,('chain',rho,sigma,crho,eps,desc))
    return tuple(res.items())
def node_items(rho,depth): return node(rho,depth)
if __name__=='__main__':
    maxd=int(sys.argv[1]) if len(sys.argv)>1 else 3
    # hub: three branches at roles a,b,c of mu. branch = child with role crho toward mu, attached at mu-role r.
    def branch(rmu,depth):
        out={}
        for crho in 'abc':
            for (TLc,TUc),desc in node(crho,depth):
                for eps in itertools.product((0,1),repeat=4):
                    R,LP,LC,E=cfgs[(rmu,crho,eps)]
                    if desc[0]=='leaf':
                        e=desc[2]
                        if any(E[(p,e)]==1 for p in [x for x in 'abc' if x!=rmu]): continue
                    CL=LP & img(R,TLc&LC); CU=LP & img(R,TUc&LC)
                    out.setdefault((CL,CU),(rmu,crho,eps,desc))
        return out
    B={r:branch(r,maxd) for r in 'abc'}
    for r in 'abc': print('branch role',r,'distinct (locked,unlocked) images:',len(B[r]))
    found=[]
    for (a1,a2) in B['a']:
        for (b1,b2) in B['b']:
            for (c1,c2) in B['c']:
                L=[a1,b1,c1]; U=[a2,b2,c2]
                if not (a2&b2&c2): continue     # tight: unlocked feasible
                ok=True
                for i in range(3):
                    s=set(D)
                    for j in range(3):
                        if j!=i: s&=L[j]
                    s&=U[i]
                    # pair locked (j!=i), i unlocked
                    if not s: ok=False
                # pairs: each pair locked, third unlocked
                if not ok: continue
                if a1&b1&c1: continue
                found.append(((a1,a2),(b1,b2),(c1,c2)))
    print('hub 3-cores found:',len(found))
    f=lambda s:''.join(sorted(s))
    for x in found[:5]: print([(f(p),f(q)) for p,q in x])
