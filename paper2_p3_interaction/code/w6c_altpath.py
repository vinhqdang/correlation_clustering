"""Path version of the alternating chain: relation between 'end slot w_0 absorbed into K of unit 1' and 'end slot w_m absorbed into K of unit m'
for all-locked / one-released anchors."""
import sys, itertools, numpy as np
from w6c_lib import *
from w6c_altchain import build
def rel(units):
    A, H, xs = build(units, ring=False)
    # recover vertex ids: rebuild index by re-running build internals
    import w6c_altchain as ac
    idx = {}
    # replicate naming to find w_0, w_m and K0 of first/last unit
    U = len(units)
    # build() creates vertices in order; find them by re-building with a hook
    return A, H, xs
def names(units):
    # reproduce idx mapping
    idx = {}
    def V(name):
        if name not in idx: idx[name] = len(idx)
        return idx[name]
    for i, t in enumerate(units):
        for v in (0, 1, 4, 5): V(('K', i, v))
        V(('w', i)); V(('w', i + 1))
        if t == 'O': V(('K', i, 7))
    return idx
def test(units):
    A, H, xs = build(units, ring=False); idx = names(units)
    opt = ilp_opt(A)[0]
    if opt != len(H): return ('notight', opt - len(H))
    U = len(units); N = len(xs)
    w0 = idx[('w', 0)]; wm = idx[('w', U)]; k0 = idx[('K', 0, 0)]; km = idx[('K', U - 1, 0)]
    def q(abs0, absm, locked):
        fx = {}
        fx[key(w0, k0)] = int(abs0); fx[key(wm, km)] = int(absm)
        for x in locked:
            for t in range(len(A)):
                if t != x: fx[key(x, t)] = 0
        r = ilp_opt(A, fixed=fx)
        return r is not None and r[0] == len(H)
    def R(locked): return {(a, b) for a in (0, 1) for b in (0, 1) if q(a, b, locked)}
    allL = R(xs); rels = [R([x for x in xs if x != xs[j]]) for j in range(N)]
    return allL, rels
for units in sys.argv[1:]:
    r = test(units)
    print(units, r if r[0] == 'notight' else ('allLocked', sorted(r[0]), 'oneReleased', [sorted(x) for x in r[1]]), flush=True)
