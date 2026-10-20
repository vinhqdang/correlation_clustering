"""w11_maj.py -- finite classification of the constraint scopes of Sigma_0 for c<=2 and closure under conservative majorities.
Scopes: triangle scope (A,B,C pairwise sharing one vertex) and edge scope (A,B sharing one vertex).
A scope relation = conjunction of all triples (inside the scope's vertices) that contain no pair owned by an outside member
(an outside-owned pair is marked X; its triples belong to another scope); remote pairs cannot occur inside a scope.
States 0,1,2 = L (edit ab), R (edit bc), M (edit ac) for member (a,b,c), centre b."""
import itertools, sys
from collections import Counter

# togetherness after edit: state -> (ab, bc, ac)
TOG = {0: (0, 1, 0), 1: (1, 0, 0), 2: (1, 1, 1)}

def member_pairs(m):
    a, b, c = m
    return [((a, b), 0), ((b, c), 1), ((a, c), 2)]

def relation(members, hard, nstates=3):
    """members: list of (a,b,c) tuples of vertex ids (centre b). hard: dict frozenset(pair)->0/1/'X'. Returns frozenset of state tuples."""
    verts = sorted(set(v for m in members for v in m))
    soft = {}
    for i, m in enumerate(members):
        for p, k in member_pairs(m):
            soft[frozenset(p)] = (i, k)
    trip = []
    for t in itertools.combinations(verts, 3):
        ps = [frozenset(q) for q in itertools.combinations(t, 2)]
        if any(hard.get(p) == 'X' for p in ps): continue
        trip.append(ps)
    res = set()
    for s in itertools.product((0, 1, 2), repeat=len(members)):
        ok = True
        for ps in trip:
            cnt = 0
            for p in ps:
                if p in soft:
                    i, k = soft[p]
                    cnt += TOG[s[i]][k]
                else:
                    cnt += hard[p]
            if cnt == 2:
                ok = False; break
        if ok: res.add(s)
    return frozenset(res)

def centre_variants(vs):
    """member on vertex list vs: 3 choices of centre -> (a,b,c)"""
    return [(vs[(i + 1) % 3], vs[i], vs[(i + 2) % 3]) for i in range(3)]

def all_relations():
    rels = {}
    # triangle scope. vertices: sAB=0, sAC=1, sBC=2, pa=3, pb=4, pc=5
    A = [0, 1, 3]; B = [0, 2, 4]; C = [1, 2, 5]
    cross_sp = [(0, 5), (1, 4), (2, 3)]   # shared-private non-member pairs
    priv = [(3, 4), (3, 5), (4, 5)]
    for ma in centre_variants(A):
        for mb in centre_variants(B):
            for mc in centre_variants(C):
                for bits in itertools.product((0, 1), repeat=3):
                    for pv in itertools.product((0, 1, 'X'), repeat=3):
                        hard = {}
                        for p, b in zip(cross_sp, bits): hard[frozenset(p)] = b
                        for p, b in zip(priv, pv): hard[frozenset(p)] = b
                        R = relation([ma, mb, mc], hard)
                        rels.setdefault(R, []).append(('tri', ma, mb, mc, bits, pv))
    # edge scope: A=(q=0,1,2), B=(q=0,3,4)
    cross = [(1, 3), (1, 4), (2, 3), (2, 4)]
    for ma in centre_variants([0, 1, 2]):
        for mb in centre_variants([0, 3, 4]):
            for cv in itertools.product((0, 1, 'X'), repeat=4):
                hard = {frozenset(p): b for p, b in zip(cross, cv)}
                R = relation([ma, mb], hard)
                rels.setdefault(R, []).append(('edge', ma, mb, cv))
    return rels

def maj_ops():
    """all 729 uniform majority operations on {0,1,2}: free values on the 6 all-distinct triples"""
    perms = list(itertools.permutations((0, 1, 2)))
    for vals in itertools.product((0, 1, 2), repeat=6):
        tab = {}
        for x, y, z in itertools.product((0, 1, 2), repeat=3):
            if x == y or x == z: tab[(x, y, z)] = x
            elif y == z: tab[(x, y, z)] = y
        for p, v in zip(perms, vals): tab[p] = v
        yield vals, tab

def closed(R, tab):
    Rl = list(R)
    for u in Rl:
        for v in Rl:
            for w in Rl:
                if tuple(tab[(u[i], v[i], w[i])] for i in range(len(u))) not in R:
                    return False
    return True

if __name__ == '__main__':
    rels = all_relations()
    print('distinct relations:', len(rels))
    byscope = Counter()
    for R, cf in rels.items():
        for c in cf: byscope[c[0]] += 0
    surv = []
    ops = list(maj_ops())
    # M-majority = all-distinct -> 2
    Mvals = tuple(2 for _ in range(6))
    Mtab = dict(ops)[0] if False else None
    for vals, tab in ops:
        if vals == Mvals: Mtab = tab
    bad = [(R, cf[0]) for R, cf in rels.items() if not closed(R, Mtab)]
    print('relations not closed under M-majority:', len(bad))
    for R, cf in bad[:10]:
        print('  ', cf, 'size', len(R))
    # which configurations fail
    kinds = Counter()
    for R, cf in bad:
        for c in rels[R]:
            if c[0] == 'tri': kinds[('tri', sum(1 for x in c[5] if x == 'X'))] += 1
            else: kinds[('edge', sum(1 for x in c[3] if x == 'X'))] += 1
    print('failing configurations by (scope, #X):', dict(kinds))
    # per-#X totals
    tot = Counter()
    for R, cfs in rels.items():
        for c in cfs:
            if c[0] == 'tri': tot[('tri', sum(1 for x in c[5] if x == 'X'))] += 1
            else: tot[('edge', sum(1 for x in c[3] if x == 'X'))] += 1
    print('total configurations by (scope,#X):', dict(tot))
    # majorities preserving all relations without X=3 (K4) configurations
    def relevant(cf):
        return not (cf[0] == 'tri' and sum(1 for x in cf[5] if x == 'X') == 3)
    Rs_noK4 = [R for R, cfs in rels.items() if any(relevant(c) for c in cfs)]
    Rs_all = list(rels)
    for name, Rs in (('all scopes incl. K4-context', Rs_all), ('without K4-context triangles', Rs_noK4)):
        good = [vals for vals, tab in ops if all(closed(R, tab) for R in Rs)]
        print(name, '#relations', len(Rs), 'majorities preserving all:', len(good), good[:5])
