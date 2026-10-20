import itertools
from a12_lib import *
S5=['L','R','M','O','X']
def allowed(Sadj, triples=(0,1,2)):
    """member vertices a,b,c=0,1,2 ; outside t=3 adjacent to Sadj; triples: 0={a,b,t},1={b,c,t},2={a,c,t}"""
    out=[]
    T=[(0,1,3),(1,2,3),(0,2,3)]
    for s in S5:
        E=apply_state({(0,1),(1,2)},(0,1,2),ST[s])|{(x,3) for x in Sadj}
        E={K(*e) for e in E}
        if all(not isP3(E,T[i]) for i in triples): out.append(s)
    return frozenset(out)
def isP3(E,T):
    a,b,c=T; return (K(a,b) in E)+(K(b,c) in E)+(K(a,c) in E)==2
# table of full-group lists
for k in range(4):
    for Sadj in itertools.combinations(range(3),k):
        print(''.join('abc'[i] for i in Sadj) or '{}',sorted(allowed(Sadj)))
# claim 7 over all single-triple-subset lists and intersections of any number of them
lists=set()
for Sadj in [S for k in range(4) for S in itertools.combinations(range(3),k)]:
    for k in range(1,4):
        for tr in itertools.combinations(range(3),k): lists.add(allowed(Sadj,tr))
print(len(lists),'distinct unary lists from one outside vertex (any subset of its 3 triples)')
viol=0
fam=list(lists)
for r in (1,2,3):
    for comb in itertools.combinations(fam,r):
        I=frozenset(S5)
        for l in comb: I&=l
        if 'O' not in I and 'L' in I and 'R' in I and (len(comb)==1 or True):
            viol+=1
print('lists/intersections (<=3 vertices) excluding O but keeping both L and R:',viol)
