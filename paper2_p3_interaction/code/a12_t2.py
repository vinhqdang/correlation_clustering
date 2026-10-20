"""Claim (2): relations between members on ZE domain {L,R,M}: closure under majority operations; also on 5-state domain."""
import itertools
from a12_lib import *
ST3=['L','R','M']; ST5=['L','R','M','O','X']
def yv(s,u,v):
    """final status of pair (u,v) positions in {0,1,2}=a,b,c : pairs ab=0, bc=1, ac=2"""
    y=ST[s]; pair={frozenset((0,1)):0,frozenset((1,2)):1,frozenset((0,2)):2}[frozenset((u,v))]
    return y[pair]
def rel(qa,pa,qb,pb,eps,dom):
    R=set()
    for sa in dom:
        for sb in dom:
            e1=yv(sa,pa,qa); e2=yv(sb,qb,pb); e3=eps
            if e1+e2+e3==2: continue   # induced P3 on {p,q,r}
            R.add((sa,sb))
    return frozenset(R)
def all_rels(dom):
    out=set()
    for qa in range(3):
      for pa in range(3):
        if pa==qa: continue
        for qb in range(3):
          for pb in range(3):
            if pb==qb: continue
            for eps in (0,1): out.add(rel(qa,pa,qb,pb,eps,dom))
    return out
def closed(R,m):
    for t1,t2,t3 in itertools.product(R,repeat=3):
        if (m(t1[0],t2[0],t3[0]),m(t1[1],t2[1],t3[1])) not in R: return False
    return True
def majorities(dom):
    trip=[t for t in itertools.permutations(dom,3)] if len(dom)==3 else None
    # majority ops: m(x,x,y)=m(x,y,x)=m(y,x,x)=x; free on all-distinct triples
    dist=[t for t in itertools.product(dom,repeat=3) if len(set(t))==3]
    for choice in itertools.product(dom,repeat=len(dist)):
        tab=dict(zip(dist,choice))
        def m(x,y,z,tab=tab):
            if x==y or x==z: return x
            if y==z: return y
            return tab[(x,y,z)]
        yield m
for dom in (ST3,):
    Rs=all_rels(dom)
    print('domain',dom,'distinct relations',len(Rs))
    ms=list(majorities(dom)); print('majorities',len(ms))
    bad=0; badR=set()
    for R in Rs:
        for m in ms:
            if not closed(R,m): bad+=1; badR.add(R)
    print('violations (rel,maj) pairs:',bad,'rels violating:',len(badR))
# extended domain: conservative majorities (median wrt a linear order) closure
dom=ST5
Rs=all_rels(dom); print('5-state relations',len(Rs))
for order in itertools.permutations(dom):
    pos={s:i for i,s in enumerate(order)}
    def med(x,y,z): return sorted((x,y,z),key=lambda s:pos[s])[1]
    if all(closed(R,med) for R in Rs): print('5-state median closed for order',order); break
else: print('no median order closes all 5-state relations')
# is there ANY majority on 5-state closing everything? (check only 'pointwise per-coordinate' not feasible; skip)
print('--- listing')
for qa in range(3):
  for pa in range(3):
    if pa==qa: continue
    for qb in range(3):
      for pb in range(3):
        if pb==qb: continue
        for eps in (0,1):
            R=rel(qa,pa,qb,pb,eps,ST3)
            if qa==qb and pa<=min(set(range(3))-{qa}) and pb<=min(set(range(3))-{qb}) or True:
                pass
seen=set()
for qa in range(3):
  for pa in range(3):
    if pa==qa: continue
    for qb in range(3):
      for pb in range(3):
        if pb==qb: continue
        for eps in (0,1):
            R=rel(qa,pa,qb,pb,eps,ST3)
            if R in seen: continue
            seen.add(R)
            ms=majorities(ST3)
            nb=sum(1 for m in ms if not closed(R,m))
            print((qa,pa,qb,pb,eps),sorted(R),'bad majorities',nb)
