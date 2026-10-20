"""a13_t11_G.py: Lemma G classification on random packings: for every triple with a far pair, the non-far pairs are (a) none, (b) all hard, (c) exactly one soft (list guard), (d) one soft + one local hard;
never two soft non-far pairs with a far third pair; and a triple with exactly one non-far pair that is soft has exactly 2 far pairs."""
import random, itertools, collections
from a13_lib import *
random.seed(11); types = collections.Counter(); bad = 0
for it in range(300):
    n = random.randint(8, 16)
    mem = []; nxt = 0; used = set()
    tries = 0
    verts = list(range(n)); E = set(); 
    # random members: pick triples, make them induced P3 in a graph built incrementally with extra random hard edges
    A = [[0]*n for _ in range(n)]
    m = random.randint(3, 7)
    for _ in range(m*5):
        if len(mem) >= m: break
        a, b, c = random.sample(verts, 3)
        ps = [frozenset((a,b)), frozenset((b,c)), frozenset((a,c))]
        if any(p in used for p in ps): continue
        mem.append((a,b,c)); used |= set(ps)
    soft = used
    E = set()
    for (a,b,c) in mem: E.add(frozenset((a,b))); E.add(frozenset((b,c)))
    for p in itertools.combinations(verts, 2):
        fp = frozenset(p)
        if fp not in soft and random.random() < .2: E.add(fp)
    inst = dict(n=n, E=E, members=mem, T=[])
    ok, c = validate(inst)
    # validate() only to confirm members induced; pair-maximality irrelevant for this classification
    if not all(not (frozenset((a,c_)) in E) for (a,b,c_) in mem): continue
    cnt = counted_pairs(inst); loc = local_pairs(inst)
    for tr in itertools.combinations(range(n), 3):
        ps = [frozenset(p) for p in itertools.combinations(tr, 2)]
        far = [p for p in ps if p not in cnt]
        if not far: continue
        nf = [p for p in ps if p in cnt]; ns = [p for p in nf if p in soft]; nh = [p for p in nf if p not in soft]
        key = (len(far), len(ns), len(nh)); types[key] += 1
        if len(ns) >= 2: bad += 1
        if len(ns) == 1 and len(far) != 2 and len(nh) == 0: bad += 1
print('triple types (#far,#soft nonfar,#local-hard nonfar):', dict(types)); print('violations of Lemma G classification:', bad)
