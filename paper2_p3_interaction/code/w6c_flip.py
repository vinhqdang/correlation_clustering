"""Search 3-member chains M1 - M2 - M3 (7 vertices) for a 'flip step':  M2 has an anchor vertex x (not shared, no hard edge),
locked state a (x isolated: L if x is the first end, R if last end).  Want: relation restricted to mid=a is a bijection D1->D3
(pass), and for some other mid state b it is a different bijection between the same 2-sets (flip)."""
import sys, itertools, numpy as np, pickle
from w6c_rel import *
names = 'LRM'
def build(p1, q1, q2, p3):
    # M1 = (0,1,2); M2 uses vertices: shared v=p1-th vertex of M1, shared w with M3
    M1 = [0, 1, 2]
    M2 = [None] * 3; M2[q1] = M1[p1]; M2[q2] = 3
    x = 4; M2[3 - q1 - q2] = x
    M3 = [None] * 3; M3[p3] = 3; others = [5, 6]
    k = 0
    for i in range(3):
        if M3[i] is None: M3[i] = others[k]; k += 1
    return [tuple(M1), tuple(M2), tuple(M3)], x, 3 - q1 - q2
found = {}
for p1 in range(3):
  for q1 in range(3):
    for q2 in range(3):
      if q2 == q1: continue
      for p3 in range(3):
        members, x, xpos = build(p1, q1, q2, p3)
        if xpos == 1: continue
        res, hard = relations(7, members)
        lock = 0 if xpos == 0 else 1           # L (del ab) isolates first end; R isolates last end
        for rel, pats in res.items():
            # patterns with no hard edge at x: need check on pattern bits; approximate: filter below
            Ra = {(s[0], s[2]) for s in rel if s[1] == lock}
            for b in range(3):
                if b == lock: continue
                Rb = {(s[0], s[2]) for s in rel if s[1] == b}
                if len(Ra) == 2 and len(Rb) == 2 and len({a for a, _ in Ra}) == 2 and len({c for _, c in Ra}) == 2 \
                   and {a for a, _ in Ra} == {a for a, _ in Rb} and {c for _, c in Ra} == {c for _, c in Rb} \
                   and len({a for a, _ in Rb}) == 2 and len({c for _, c in Rb}) == 2 and Ra != Rb:
                    # check some pattern realising rel has no hard edge at x
                    hidx = {p: i for i, p in enumerate(hard)}
                    good = []
                    for j in pats:
                        bad = any((j >> hidx[key(x, t)]) & 1 for t in range(7) if t != x and key(x, t) in hidx)
                        if not bad: good.append(j)
                    if good:
                        found.setdefault((p1, q1, q2, p3), []).append((sorted(Ra), b, sorted(Rb), good[:2]))
print(len(found))
for k, v in found.items(): print(k, v[:2])
pickle.dump(found, open(outpath('flip_found.pkl'), 'wb'))
