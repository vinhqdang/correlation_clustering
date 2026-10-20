"""w9b_sigma.py -- vectorised ZE-set Sigma_0 of an instance (far pairs dropped) and list-core numbers."""
import itertools
import numpy as np
from w9b_lib import Inst, key

XT = np.array([[0, 1, 1], [1, 0, 1], [0, 0, 1]])   # XT[role ab/bc/ac][state L,R,M]  -> togetherness

def sigma(inst, drop_far=True, extra_lists=None):
    """returns array (K, m) of ZE state vectors, constraints = all triples without a dropped far pair"""
    n = inst.n; m = len(inst.H)
    if m == 0: return np.zeros((1, 0), int)
    S = np.array(list(itertools.product((0, 1, 2), repeat=m)), dtype=np.int8) if m <= 11 else None
    ok = np.ones(len(S), bool)
    far = set(inst.far_pairs()) if drop_far else set()
    for x, y, z in itertools.combinations(range(n), 3):
        ps = [key(x, y), key(y, z), key(x, z)]
        if any(p in far for p in ps): continue
        tot = np.zeros(len(S), np.int8)
        anysoft = False
        for p in ps:
            o = inst.owner.get(p)
            if o is None:
                tot += 1 if inst.adj[p[0]][p[1]] else 0
            else:
                anysoft = True
                tot += XT[o[1]][S[:, o[0]]]
        if not anysoft:
            if tot.max() == 2 and tot.min() == 2: ok[:] = False
            continue
        ok &= (tot != 2)
    return S[ok]

LISTS = [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2)]

def list_masks(Sg):
    """for each (member, list) the bitmask (python int) over rows of Sg satisfying it; drop trivial/duplicate"""
    K, m = Sg.shape
    items = []; seen = {}
    full = (1 << K) - 1
    for i in range(m):
        for L in LISTS:
            sat = np.isin(Sg[:, i], L)
            mask = 0
            idx = np.nonzero(sat)[0]
            for j in idx: mask |= 1 << int(j)
            if mask == full: continue
            items.append((i, L, mask))
    return items, full

def max_core(Sg, maxk=8, singleton_only=False, return_family=False):
    """max size of a minimal inconsistent family of lists (one list per member in the family).  DFS with strictly decreasing AND."""
    items, full = list_masks(Sg)
    if singleton_only: items = [it for it in items if len(it[1]) == 1]
    best = [0, None]
    # sort so that DFS is canonical
    N = len(items)
    def rec(start, cur, fam, used):
        for j in range(start, N):
            i, L, mk = items[j]
            if i in used: continue
            nxt = cur & mk
            if nxt == cur: continue            # not essential
            if nxt == 0:
                # candidate inconsistent family; check minimality
                f2 = fam + [j]
                if len(f2) > best[0]:
                    good = True
                    for t in range(len(f2)):
                        a = full
                        for u in range(len(f2)):
                            if u != t: a &= items[f2[u]][2]
                        if a == 0: good = False; break
                    if good:
                        best[0] = len(f2); best[1] = [(items[u][0], items[u][1]) for u in f2]
                continue
            if len(fam) + 1 < maxk:
                rec(j + 1, nxt, fam + [j], used | {i})
    rec(0, full, [], frozenset())
    return (best[0], best[1]) if return_family else best[0]
