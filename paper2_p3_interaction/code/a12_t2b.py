import itertools
from a12_lib import *
from a12_t2 import yv, majorities, closed, ST3, ST5
def relc(qa,qb,bits,dom):
    pas=[x for x in range(3) if x!=qa]; pbs=[x for x in range(3) if x!=qb]
    R=set()
    for sa in dom:
        for sb in dom:
            ok=True
            for i,pa in enumerate(pas):
                for j,pb in enumerate(pbs):
                    if yv(sa,pa,qa)+yv(sb,qb,pb)+bits[2*i+j]==2: ok=False
            if ok: R.add((sa,sb))
    return frozenset(R)
for dom in (ST3,ST5):
    Rs=set(relc(qa,qb,b,dom) for qa in range(3) for qb in range(3) for b in itertools.product((0,1),repeat=4))
    print(dom,'distinct conjunction relations',len(Rs),'sizes',sorted(len(r) for r in Rs)[:0] or sorted(set(len(r) for r in Rs)))
    if len(dom)==3:
        ms=list(majorities(dom)); bad=sum(1 for R in Rs for m in ms if not closed(R,m)); print(' violations over all 729 majorities:',bad)
    else:
        # test median closure for all orders, and a search for any uniform majority is infeasible (5^60); report medians
        found=0
        for order in itertools.permutations(dom):
            pos={s:i for i,s in enumerate(order)}
            med=lambda x,y,z:sorted((x,y,z),key=lambda s:pos[s])[1]
            if all(closed(R,med) for R in Rs): found+=1
        print(' median orders closing all 5-state conj. relations:',found,'of 120')
        # which relations fail some median / the notes claim extended domain not majority closed
        # Boolean-majority style: all-distinct triples -> try per-triple constraint by exhaustive over relations using CSP of choices
