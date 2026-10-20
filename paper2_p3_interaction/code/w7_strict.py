"""Relation-level lemma: for every edge configuration (R, LamP, LamC) and all X <= Y subsets of D:
X' = LamP & img(R, X & LamC), Y' similarly; if X' != Y' then |X'| <= 1."""
import itertools, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from w7_hub3 import cfgs, img, D, S
allsub=[frozenset(c) for r in range(0,4) for c in itertools.combinations(S,r)]
viol=[];tested=0;strict=0
for key,(R,LP,LC,E) in cfgs.items():
    for X in allsub:
        for Y in allsub:
            if not X<=Y: continue
            X2=LP&img(R,X&LC); Y2=LP&img(R,Y&LC); tested+=1
            if X2!=Y2:
                strict+=1
                if len(X2)>=2: viol.append((key,sorted(X),sorted(Y),sorted(X2),sorted(Y2)))
print('configs',len(cfgs),'tested (X<=Y)',tested,'strict',strict,'violations of |X\'|<=1:',len(viol))
print(viol[:3])
