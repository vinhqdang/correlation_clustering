"""Wave 10: the reduction Grid Tiling (<=)  ->  CEaMP, builder (thin guest model or thick), semantics tests.

Grid Tiling with <=: k, n, S[i][j] subset of [n]^2 (1-based values).  Find s_ij=(a_ij,b_ij) in S_ij with a_ij <= a_{i,j+1}, b_ij <= b_{i+1,j}.
Chains: A_ij (tiles in lexicographic order (a,b): block a, slot b), B_ij (order (b,a): block b, slot a).
Chain layout (all chains identical): groups u = 0..U-1, U = 1 + n(n+1); group u has positions (2u, 2u+1), host s_{4u+2};
   group 0 and every separator group: end/separator pad; slot groups (alpha-1)(n+1)+1+(slot-1) for block alpha=1..n; separator after block alpha at group alpha(n+1).
Slot group of a tile in S_ij: lock guest + anchor (tunnel to the partner chain); slot group of a non-tile: permanent forbid (pad).
Comparators: separator alpha in 1..n-1 between consecutive chains of a row (A) / column (B).
"""
import sys, os, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w10_lib import *


def layout(n):
    U = 1 + n * (n + 1)
    m = 2 * U
    slot_group = {}
    sep_group = {}
    for alpha in range(1, n + 1):
        for s in range(1, n + 1):
            slot_group[(alpha, s)] = (alpha - 1) * (n + 1) + s
        sep_group[alpha] = alpha * (n + 1)
    return U, m, slot_group, sep_group


def build_gt(k, n, S, p, thin_mode=True, tunnels=True):
    U, m, slot_group, sep_group = layout(n)
    W = World(p, thin_mode=thin_mode)
    cells = [(i, j) for i in range(k) for j in range(k)]
    chains = {}
    for (i, j) in cells:
        chains[('A', i, j)] = W.add_chain(('A', i, j), m)
        chains[('B', i, j)] = W.add_chain(('B', i, j), m)
    zname = lambda tag, u: ('z', tag, u)
    meta = dict(U=U, m=m, slot_group=slot_group, sep_group=sep_group, anchors={}, tunnels=[], k=k, n=n, S=S, chains=chains, locks={})
    # which groups have a lock guest
    tile_groups = {tag: {} for tag in chains}
    for (i, j) in cells:
        for (a, b) in S[i][j]:
            tile_groups[('A', i, j)][slot_group[(a, b)]] = (a, b)
            tile_groups[('B', i, j)][slot_group[(b, a)]] = (a, b)
    sep_set = set(sep_group.values())
    for tag, ch in chains.items():
        for u in range(U):
            host = ('s', 4 * u + 2)
            F = {2 * u, 2 * u + 1}
            if u in tile_groups[tag]:
                # permanent guest allows everything, lock guest forbids F when locked
                ch.perm_guest(zname(tag, u), host, set())
                x, y = ('x', tag, u), ('y', tag, u)
                lock_guest(ch, y, x, zname(tag, u), host, F)
                meta['anchors'][(tag, u)] = x
            else:
                eff = ch.perm_guest(zname(tag, u), host, F)
                assert eff == F, (tag, u, eff)
    # tunnels
    for (i, j) in cells:
        for (a, b) in S[i][j]:
            xa = ('x', ('A', i, j), slot_group[(a, b)])
            xb = ('x', ('B', i, j), slot_group[(b, a)])
            meta['tunnels'].append(((i, j), (a, b), xa, xb))
            if tunnels:
                W.B.hard(xa, xb)
    # comparators
    cmpl = []
    for i in range(k):
        for j in range(k - 1):
            for alpha in range(1, n):
                cmpl.append((('A', i, j), ('A', i, j + 1), alpha))
    for j in range(k):
        for i in range(k - 1):
            for alpha in range(1, n):
                cmpl.append((('B', i, j), ('B', i + 1, j), alpha))
    for (c1, c2, alpha) in cmpl:
        u = sep_group[alpha]
        host = ('s', 4 * u + 2)
        F1 = {a for a in range(m) if a >= 2 * u}
        F2 = {a for a in range(m) if a <= 2 * u + 1}
        y = ('c', c1, c2, alpha)
        effA, effB = cmp_guest(chains[c1], host, F1, chains[c2], host, F2, y, zname(c1, u), zname(c2, u))
        assert effA == {a for a in range(m) if a >= 2 * u} and effB == {a for a in range(m) if a <= 2 * u + 1}, (effA, effB)
    meta['cmp'] = cmpl
    return W, meta


def gt_solvable(k, n, S):
    cells = [(i, j) for i in range(k) for j in range(k)]
    sols = []
    for choice in itertools.product(*[sorted(S[i][j]) for (i, j) in cells]):
        t = {c: e for c, e in zip(cells, choice)}
        ok = True
        for i in range(k):
            for j in range(k):
                if j + 1 < k and t[(i, j)][0] > t[(i, j + 1)][0]:
                    ok = False
                if i + 1 < k and t[(i, j)][1] > t[(i + 1, j)][1]:
                    ok = False
        if ok:
            sols.append(t)
    return sols


def cut_choice(meta, t):
    """tunnel edges to delete for tile choice t (dict cell->tile)"""
    drops = []
    for (cell, e, xa, xb) in meta['tunnels']:
        if t[cell] == e:
            drops.append((xa, xb))
    return drops


if __name__ == '__main__':
    import random
    k, n = int(sys.argv[1]), int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    dens = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
    rng = random.Random(seed)
    S = [[set() for _ in range(k)] for _ in range(k)]
    for i in range(k):
        for j in range(k):
            for a in range(1, n + 1):
                for b in range(1, n + 1):
                    if rng.random() < dens:
                        S[i][j].add((a, b))
            if not S[i][j]:
                S[i][j].add((rng.randint(1, n), rng.randint(1, n)))
    print('S', S)
    t0 = time.time()
    W, meta = build_gt(k, n, S, 2, thin_mode=True)
    I = ThinInst(W)
    print('thin n', I.n, '|H|', len(I.members), I.valid, round(time.time() - t0, 1), flush=True)
    sols = gt_solvable(k, n, S)
    print('GT solutions', len(sols))
    cells = [(i, j) for i in range(k) for j in range(k)]
    agree = 0; tot = 0
    for choice in itertools.product(*[sorted(S[i][j]) for (i, j) in cells]):
        t = {c: e for c, e in zip(cells, choice)}
        drops = cut_choice(meta, t)
        lab = I.solve(drop_hard=drops)
        exp = t in sols
        tot += 1
        agree += ((lab is not None) == exp)
        if (lab is not None) != exp:
            print('MISMATCH', t, lab is not None, exp, flush=True)
    print('agree', agree, 'of', tot, round(time.time() - t0, 1))


def layout2(n, extra=0):
    U, m, slot_group, sep_group = layout(n)
    if extra:
        U += extra; m = 2 * U
    return U, m, slot_group, sep_group


def build_general(n, p, tiles, cmps, tunnels, thin_mode=False, extra=0):
    """tiles: dict tag -> set of slot-group indices carrying a lock guest+anchor; every other group is a permanent pad.
    cmps: list of (tag1, tag2, alpha); tunnels: list of ((tag,u),(tag',u')) or ((tag,u), None) = private stub.
    All chains share the layout of Grid Tiling with parameter n (U=1+n(n+1)+extra groups)."""
    U, m, slot_group, sep_group = layout2(n, extra)
    W = World(p, thin_mode=thin_mode)
    chains = {tag: W.add_chain(tag, m) for tag in tiles}
    zname = lambda tag, u: ('z', tag, u)
    for tag, ch in chains.items():
        for u in range(U):
            host = ('s', 4 * u + 2)
            F = {2 * u, 2 * u + 1}
            if u in tiles[tag]:
                ch.perm_guest(zname(tag, u), host, set())
                lock_guest(ch, ('y', tag, u), ('x', tag, u), zname(tag, u), host, F)
            else:
                eff = ch.perm_guest(zname(tag, u), host, F)
                assert eff == F
    for (a, b) in tunnels:
        xa = ('x',) + a
        if b is None:
            W.B.hard(xa, ('stub',) + a)
        else:
            W.B.hard(xa, ('x',) + b)
    for (c1, c2, alpha) in cmps:
        u = sep_group[alpha]
        host = ('s', 4 * u + 2)
        F1 = {a for a in range(m) if a >= 2 * u}
        F2 = {a for a in range(m) if a <= 2 * u + 1}
        cmp_guest(chains[c1], host, F1, chains[c2], host, F2, ('c', c1, c2, alpha), zname(c1, u), zname(c2, u))
    meta = dict(U=U, m=m, slot_group=slot_group, sep_group=sep_group, chains=chains)
    return W, meta
