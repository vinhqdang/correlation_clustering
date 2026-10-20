"""a13_t6_maj.py: (i) Theorem A at graph level: for two members sharing exactly one vertex (all orientations, all 16 cross-pair statuses),
the full conjunction of triple constraints on their 5 vertices is preserved by ALL 729 uniform majorities (all-distinct triples -> arbitrary value,
same function in every coordinate).  (ii) Lemma H test: random majority-closed relations of arity 4,5 over a 3-element domain, all list systems."""
import itertools, random
def pres(R, arity=None):
    """True iff every uniform majority preserves R (set of tuples over {0,1,2}); returns offending triple if not"""
    R = list(R); Rs = set(R); k = len(R[0]) if R else 0
    for x, y, z in itertools.product(R, repeat=3):
        fixed = {}; free = []
        out = [None]*k
        for c in range(k):
            a, b, d = x[c], y[c], z[c]
            if a == b or a == d: out[c] = a
            elif b == d: out[c] = b
            else: free.append(c)
        keys = sorted(set((x[c], y[c], z[c]) for c in free))
        for vals in itertools.product(range(3), repeat=len(keys)):
            m = dict(zip(keys, vals)); o = list(out)
            for c in free: o[c] = m[(x[c], y[c], z[c])]
            if tuple(o) not in Rs: return (x, y, z, m)
    return None
def closure(S):
    S = set(S)
    while True:
        new = set()
        R = list(S)
        for x, y, z in itertools.product(R, repeat=3):
            k = len(x); out = [None]*k; free = []
            for c in range(k):
                a, b, d = x[c], y[c], z[c]
                if a == b or a == d: out[c] = a
                elif b == d: out[c] = b
                else: free.append(c)
            keys = sorted(set((x[c], y[c], z[c]) for c in free))
            for vals in itertools.product(range(3), repeat=len(keys)):
                m = dict(zip(keys, vals)); o = list(out)
                for c in free: o[c] = m[(x[c], y[c], z[c])]
                if tuple(o) not in S: new.add(tuple(o))
        if not new: return S
        S |= new
# ---------- (i)
STATE = (0, 1, 2)  # L,R,M
def rel(A, B, cross):   # A=(a,b,c), B=(a,b,c) vertex labels; cross: dict pair->0/1
    verts = sorted(set(A) | set(B)); R = set()
    for sa in STATE:
        for sb in STATE:
            adj = {}
            for (a, b, c), s in ((A, sa), (B, sb)):
                adj[frozenset((a, b))] = 0 if s == 0 else 1
                adj[frozenset((b, c))] = 0 if s == 1 else 1
                adj[frozenset((a, c))] = 1 if s == 2 else 0
            ok = True
            for x, y, z in itertools.combinations(verts, 3):
                ps = [frozenset((x, y)), frozenset((y, z)), frozenset((x, z))]
                t = sum(adj[p] if p in adj else cross[p] for p in ps)
                if t == 2: ok = False; break
            if ok: R.add((sa, sb))
    return frozenset(R)
rels = {}
for A in itertools.permutations((0, 1, 2)):
    for B in itertools.permutations((0, 3, 4)):    # vertex 0 shared
        cp = [(1, 3), (1, 4), (2, 3), (2, 4)]
        for st in itertools.product((0, 1), repeat=4):
            cross = {frozenset(p): s for p, s in zip(cp, st)}
            rels.setdefault(rel(A, B, cross), (A, B, st))
print('distinct binary relations (incl. empty):', len(rels), ' empty ones:', int(frozenset() in rels))
bad = [(R, c) for R, c in rels.items() if R and pres(R)]
print('relations NOT preserved by some uniform majority:', len(bad))
# unary lists
for S in itertools.chain.from_iterable(itertools.combinations(range(3), r) for r in (1, 2, 3)):
    assert pres(set((s,) for s in S)) is None
print('lists preserved')
# ---------- (ii) Lemma H
import time
random.seed(5)
viol = 0; tested = 0; maxI = 0; nrel = 0
opts = [frozenset(c) for r in (1, 2) for c in itertools.combinations(range(3), r)] + [frozenset(range(3))]
t0 = time.time()
for ar in (4, 5):
    for it in range(150):
        if time.time() - t0 > 240: break
        pts = list(itertools.product(range(3), repeat=ar))
        S = set(random.sample(pts, random.randint(2, 4)))
        for _ in range(6):
            R = list(S); new = set()
            for x, y, z in itertools.product(R, repeat=3):
                out = [None]*ar; free = []
                for c in range(ar):
                    a, b, d = x[c], y[c], z[c]
                    if a == b or a == d: out[c] = a
                    elif b == d: out[c] = b
                    else: free.append(c)
                keys = sorted(set((x[c], y[c], z[c]) for c in free))
                for vals in itertools.product(range(3), repeat=len(keys)):
                    m = dict(zip(keys, vals)); o = list(out)
                    for c in free: o[c] = m[(x[c], y[c], z[c])]
                    new.add(tuple(o))
            if new <= S: break
            S |= new
            if len(S) > 40: break
        if len(S) > 40 or pres(S) is not None: continue
        nrel += 1
        for L in itertools.product(opts, repeat=ar):
            if any(all(p[i] in L[i] for i in range(ar)) for p in S): continue
            I = [i for i in range(ar) if len(L[i]) < 3]
            minimal = True
            for i in I:
                L2 = list(L); L2[i] = frozenset(range(3))
                if not any(all(p[j] in L2[j] for j in range(ar)) for p in S): minimal = False; break
            if minimal:
                tested += 1; maxI = max(maxI, len(I))
                if len(I) > 3: viol += 1
print('Lemma H: %d closed relations (arity 4,5), minimal inconsistent list systems=%d, max #coordinates=%d, violations(>3)=%d' % (nrel, tested, maxI, viol), flush=True)
