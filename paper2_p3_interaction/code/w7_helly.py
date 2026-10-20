"""Fixpoint computation for forests (no soft triangle): family of achievable (locked-image X, unlocked-image Y) pairs of a branch, as seen
from the member at which it attaches. Generous model: arbitrary static unary lists may be intersected anywhere.
Property S: whenever a lock really restricts (X != Y) the locked image has at most one element -> no hub 3-core."""
import itertools, io, contextlib, sys
with contextlib.redirect_stdout(io.StringIO()):
    from w7_hub3 import cfgs, img, D, S
allsub=[frozenset(c) for r in range(0,4) for c in itertools.combinations(S,r)]
mode=sys.argv[1] if len(sys.argv)>1 else 'single'
def leafpairs():
    out=set()
    if mode=='single':
        for t in 'LR': out.add((frozenset(t),D))
    else:
        for T in allsub:
            if 0<len(T)<3: out.add((T,D))
    return out
# family per role rho (role of the vertex through which the branch's root member attaches to its parent)
fam={r:set(leafpairs())|{(D,D)} for r in 'abc'}
def close_lists(F):
    out=set(F)
    for (X,Y) in F:
        for lam in allsub:
            out.add((X&lam,Y&lam))
    return out
changed=True; it=0
while changed:
    changed=False; it+=1
    new={r:set(fam[r]) for r in 'abc'}
    for rho in 'abc':
        for sigma in 'abc':
            if sigma==rho: continue
            for crho in 'abc':
                for (X,Y) in fam[crho]:
                    for eps in itertools.product((0,1),repeat=4):
                        R,LP,LC,E=cfgs[(sigma,crho,eps)]
                        X2=LP&img(R,X&LC); Y2=LP&img(R,Y&LC)
                        if Y2: new[rho].add((X2,Y2))
        new[rho]=close_lists(new[rho])
    for r in 'abc':
        if new[r]!=fam[r]: fam[r]=new[r]; changed=True
print('mode',mode,'fixpoint after',it,'rounds; |family| per role',{r:len(fam[r]) for r in fam})
# images seen at a hub through relation at hub vertex roles
f=lambda s:''.join(sorted(s)) or '0'
bad=set()
imgs={r:set() for r in 'abc'}
for rmu in 'abc':
    for crho in 'abc':
        for (X,Y) in fam[crho]:
            for eps in itertools.product((0,1),repeat=4):
                R,LP,LC,E=cfgs[(rmu,crho,eps)]
                X2=LP&img(R,X&LC); Y2=LP&img(R,Y&LC)
                if Y2: imgs[rmu].add((X2,Y2))
for r in 'abc':
    strict2=sorted({f(X) for X,Y in imgs[r] if X!=Y and len(X)>=2})
    print('hub vertex role',r,'strict locked images of size>=2:',strict2,'(total pairs',len(imgs[r]),')')
# hub 3-core search with all three branches + own list
found=0
for (a1,a2) in imgs['a']:
    for (b1,b2) in imgs['b']:
        for (c1,c2) in imgs['c']:
            if not (a2&b2&c2): continue
            L=[a1,b1,c1];U=[a2,b2,c2]
            if all((U[i]&L[j]&L[k]) for i,j,k in ((0,1,2),(1,0,2),(2,0,1))) and not (a1&b1&c1): found+=1
print('hub 3-cores (3 branches):',found)
