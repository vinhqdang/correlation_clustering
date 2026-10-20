"""Lemma S' (strictness) for the 40 ZE relations R on {L,R,M}^2 (both directions): for all X subset Y subset {L,R,M}: img_R(X) != img_R(Y) implies |img_R(X)| <= 1.
Then Helly-2 for lists on a forest (proof in wave9 note) and a brute-force confirmation."""
import itertools, sys
sys.argv=['x','1','1','7']
exec(open('/home/user/correlation_clustering/paper2_p3_interaction/code/w9a_helly.py').read().split("rng=random.Random")[0])
subsets=[frozenset(c) for k in range(4) for c in itertools.combinations(STS,k)]
viol=0; checked=0; strict=0
for R in RELS:
    for direction in (0,1):
        Rd=R if direction==0 else frozenset((t,s) for (s,t) in R)
        def img(X): return frozenset(s for s in STS if any((s,t) in Rd for t in X))
        for X in subsets:
            for Y in subsets:
                if X<=Y:
                    checked+=1
                    if img(X)!=img(Y):
                        strict+=1
                        if len(img(X))>1: viol+=1; print('VIOL',sorted(Rd),sorted(X),sorted(Y),sorted(img(X)),sorted(img(Y)))
print('relations',len(RELS),'checked',checked,'strict',strict,'violations',viol)
