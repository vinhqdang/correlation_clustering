"""w11_gen.py -- generators of pair-maximal instances: members (tree-like I(H), c<=2) + hard clusters (partition)."""
import random, itertools
from w11_lib import *

def random_members(m, rnd, tri=False):
    """m members whose intersection graph is a random tree (max degree 3, each vertex in <=2 members). returns (nv, H)"""
    # each member has 3 slots; centre choice random
    slots = [[None]*3 for _ in range(m)]
    nv = 0
    edges = []
    deg = [0]*m
    for i in range(1, m):
        # attach to random earlier member with free slot
        cand = [j for j in range(i) if deg[j] < 3]
        j = rnd.choice(cand)
        # choose free slot in j and in i
        sj = rnd.choice([s for s in range(3) if slots[j][s] is None])
        si = rnd.choice([s for s in range(3) if slots[i][s] is None])
        v = nv; nv += 1
        slots[j][sj] = v; slots[i][si] = v
        deg[j] += 1; deg[i] += 1
    for i in range(m):
        for s in range(3):
            if slots[i][s] is None:
                slots[i][s] = nv; nv += 1
    H = []
    for i in range(m):
        t = slots[i][:]
        rnd.shuffle(t)       # t[1] is the centre
        H.append(tuple(t))
    return nv, H

def make_inst(nv, H, npend, clusters):
    """clusters: list of cluster ids for each vertex (len nv+npend). hard pairs inside a cluster are edges."""
    n = nv + npend
    E = set()
    for (a, b, c) in H:
        E.add(key(a, b)); E.add(key(b, c))
    soft = set()
    for (a, b, c) in H:
        soft |= {key(a, b), key(b, c), key(a, c)}
    for u, v in itertools.combinations(range(n), 2):
        if key(u, v) in soft: continue
        if clusters[u] == clusters[v] and clusters[u] >= 0:
            E.add(key(u, v))
    return Inst(n, E, H)

def random_members_c2(m, rnd, density=0.7):
    """general c<=2 packing: random matching of member slots (members share <=1 vertex). returns (nv,H)."""
    slots = [(i, s) for i in range(m) for s in range(3)]
    rnd.shuffle(slots)
    vid = {}
    nv = 0
    shared_pairs = set()
    used = set()
    k = 0
    cand = slots[:]
    pairs = []
    while len(cand) >= 2:
        a = cand.pop()
        if rnd.random() > density:
            continue
        # find partner
        rnd.shuffle(cand)
        for idx, b in enumerate(cand):
            if b[0] != a[0] and frozenset((a[0], b[0])) not in shared_pairs:
                shared_pairs.add(frozenset((a[0], b[0])))
                pairs.append((a, b)); cand.pop(idx); break
    assign = {}
    for a, b in pairs:
        assign[a] = assign[b] = nv; nv += 1
    H = []
    for i in range(m):
        t = []
        for s in range(3):
            if (i, s) in assign: t.append(assign[(i, s)])
            else: t.append(nv); nv += 1
        rnd.shuffle(t)
        H.append(tuple(t))
    return nv, H

def random_members_cc(m, rnd, cmax=3, density=0.8, p3=0.4):
    """c<=cmax packing: random grouping of member slots into vertices of multiplicity 1..cmax; members share <=1 vertex."""
    slots = [(i, s) for i in range(m) for s in range(3)]
    rnd.shuffle(slots)
    shared = set()
    groups = []
    rest = slots[:]
    while rest:
        a = rest.pop()
        if rnd.random() > density: groups.append([a]); continue
        g = [a]
        target = 3 if (cmax >= 3 and rnd.random() < p3) else 2
        rnd.shuffle(rest)
        for b in rest[:]:
            if len(g) >= target: break
            if all(b[0] != x[0] and frozenset((b[0], x[0])) not in shared for x in g):
                g.append(b); rest.remove(b)
        for x in g:
            for y in g:
                if x[0] < y[0]: shared.add(frozenset((x[0], y[0])))
        groups.append(g)
    vid = {}
    for v, g in enumerate(groups):
        for x in g: vid[x] = v
    nv = len(groups)
    H = []
    for i in range(m):
        t = [vid[(i, s)] for s in range(3)]
        rnd.shuffle(t)
        H.append(tuple(t))
    return nv, H
