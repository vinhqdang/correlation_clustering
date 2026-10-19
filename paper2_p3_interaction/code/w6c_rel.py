"""Relation search: k members (a,b,c) (centre b) on n vertices.  For every pattern of the hard pairs (edge/non-edge),
the 'zero-excess relation' = set of tuples of member states (L=del ab, R=del bc, M=insert ac) that occur in a partition
of cost exactly k (no hard pair edited).  Used to find which local relations (wires, XOR-type steps, ...) are realisable."""
import itertools, numpy as np
from w6c_lib import *

def relations(n, members, hard_extra_fixed=None):
    lab, same, P = part_table(n); same_f = same.astype(np.float32)
    pidx = {p: i for i, p in enumerate(P)}
    soft = set(); fixed_edges = set()
    for (a, b, c) in members:
        soft |= pairs_of((a, b, c)); fixed_edges |= {key(a, b), key(b, c)}
    hard = [p for p in P if p not in soft]
    h = len(hard)
    Hm = ((np.arange(2 ** h)[:, None] >> np.arange(h)[None, :]) & 1).astype(np.float32)
    E = np.zeros((2 ** h, len(P)), np.float32)
    for p in fixed_edges: E[:, pidx[p]] = 1
    for i, p in enumerate(hard): E[:, pidx[p]] = Hm[:, i]
    cost = same_f @ (1 - 2 * E.T) + E.sum(1)[None, :]
    k = len(members)
    # state per member per partition
    def sig(m):
        a, b, c = m
        s = same[:, [pidx[key(a, b)], pidx[key(b, c)], pidx[key(a, c)]]].astype(int)
        out = np.full(len(s), -1)
        for code, pat in enumerate(((0, 1, 0), (1, 0, 0), (1, 1, 1))):
            out[(s == np.array(pat)).all(1)] = code
        return out
    S = np.stack([sig(m) for m in members], 1)          # (partitions, k)
    valid = (S >= 0).all(1)
    res = {}
    ok = cost == k
    for j in range(2 ** h):
        m = ok[:, j] & valid
        if m.any():
            rel = frozenset(map(tuple, S[m]))
            res.setdefault(rel, []).append(j)
    return res, hard
