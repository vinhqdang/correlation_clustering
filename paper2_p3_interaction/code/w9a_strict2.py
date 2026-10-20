"""Property I for the 40 ZE relations (both directions): for all Y, L1, L2 subsets of {L,R,M}: img(Y cap L1 cap L2) != img(Y) implies img(Y cap L1) != img(Y) or img(Y cap L2) != img(Y).
Together with strictness S' this gives: every inconsistency of a forest CSP over the family with arbitrary lists has a certificate of at most two list constraints (Helly-2)."""
import itertools, sys
sys.argv=['x','1','1','7']
exec(open('/home/user/correlation_clustering/paper2_p3_interaction/code/w9a_helly.py').read().split("rng=random.Random")[0])
subsets=[frozenset(c) for k in range(4) for c in itertools.combinations(STS,k)]
viol=0; tot=0
for R in RELS:
    for direction in (0,1):
        Rd=R if direction==0 else frozenset((t,s) for (s,t) in R)
        def img(X): return frozenset(s for s in STS if any((s,t) in Rd for t in X))
        for Y in subsets:
            for L1 in subsets:
                for L2 in subsets:
                    tot+=1
                    if (Y&L1&L2) and img(Y&L1&L2)!=img(Y) and img(Y&L1)==img(Y) and img(Y&L2)==img(Y):
                        viol+=1
                        if viol<=5: print('I-VIOL',sorted(Rd),sorted(Y),sorted(L1),sorted(L2),sorted(img(Y&L1&L2)),sorted(img(Y)))
print('checked',tot,'violations of Property I:',viol)
