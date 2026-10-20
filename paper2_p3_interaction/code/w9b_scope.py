"""w9b_scope.py -- scope relations built from REAL graphs (independent of the formula X[e][s] in w9b_poly).
A scope is a set of 1..3 mutually adjacent members on a small vertex set; hard pairs have status 0 (non-edge), 1 (edge) or None
('dropped': soft pair of an outside member, its triples belong to another scope).  The scope relation is the set of state vectors in {L,R,M}^k
(L: delete ab, R: delete bc, M: insert ac) for which no triple of the vertex set, not containing a dropped pair, is an induced P3 of the edited graph.
"""
import itertools
from functools import lru_cache

STATES = (0, 1, 2)   # L,R,M

def edited_adj(members, hard, nverts, states):
    adj = {}
    for (a, b, c), s in zip(members, states):
        adj[frozenset((a, b))] = 0 if s == 0 else 1     # edit ab for L
        adj[frozenset((b, c))] = 0 if s == 1 else 1     # edit bc for R
        adj[frozenset((a, c))] = 1 if s == 2 else 0     # insert ac for M
    return adj

def scope_relation(members, hard, nverts):
    """members: list of (a,b,c); hard: dict frozenset(pair)-> 0/1/None for every non-soft pair"""
    out = set()
    for states in itertools.product(STATES, repeat=len(members)):
        adj = edited_adj(members, hard, nverts, states)
        ok = True
        for x, y, z in itertools.combinations(range(nverts), 3):
            ps = [frozenset((x, y)), frozenset((y, z)), frozenset((x, z))]
            tot = 0; dropped = False
            for p in ps:
                if p in adj: tot += adj[p]
                else:
                    h = hard[p]
                    if h is None: dropped = True; break
                    tot += h
            if dropped: continue
            if tot == 2: ok = False; break
        if ok: out.add(states)
    return frozenset(out)

def role_orders(vs):
    """the 3 members on vertex set vs=(u,v,w): centre choices; (a,b,c) with centre b (a<->c symmetric)."""
    return list(itertools.permutations(vs))   # all 6 orientations (a,b,c): centre b, state labels L/R depend on orientation

def pair_scope_relations():
    # A on {q,p1,p2}, B on {q,r1,r2}; vertices 0=q,1=p1,2=p2,3=r1,4=r2
    rels = {}
    cross = [(1, 3), (1, 4), (2, 3), (2, 4)]
    for A in role_orders((0, 1, 2)):
        for B in role_orders((0, 3, 4)):
            for st in itertools.product((0, 1, None), repeat=4):
                hard = {frozenset(p): s for p, s in zip(cross, st)}
                R = scope_relation([A, B], hard, 5)
                rels.setdefault(R, (A, B, st))
    return rels

def tri_scope_relations():
    # A=(p,q,a) on {p=0,q=1,a=2}; B on {q=1,r=3,b=4}; C on {r=3,p=0,c=5}
    rels = {}
    hardpairs = [(2, 4), (2, 5), (4, 5), (2, 3), (4, 0), (5, 1)]   # ab, ac, bc, (a,r), (b,p), (c,q)
    for A in role_orders((0, 1, 2)):
        for B in role_orders((1, 3, 4)):
            for C in role_orders((3, 0, 5)):
                for st in itertools.product((0, 1, None), repeat=6):
                    hard = {frozenset(p): s for p, s in zip(hardpairs, st)}
                    R = scope_relation([A, B, C], hard, 6)
                    rels.setdefault(R, (A, B, C, st))
    return rels

if __name__ == '__main__':
    import pickle, time
    t = time.time()
    P = pair_scope_relations()
    print('pair scope relations (with dropped):', len(P), time.time() - t, flush=True)
    full = {R for R, (A, B, st) in P.items() if None not in st}
    print(' of which without dropped pairs:', len(full))
    T = tri_scope_relations()
    print('triangle scope relations:', len(T), time.time() - t, flush=True)
    pickle.dump((P, T), open('/tmp/w9b_scopes.pkl', 'wb'))
