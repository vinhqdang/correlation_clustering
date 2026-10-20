"""Thick (blown-up) instances: every vertex v of a thin instance becomes a proto-cluster K_p(v) (vertices (v,i), i<p, hard-glued),
every member (a,b,c) becomes the p^2 Latin members ((a,i),(b,j),(c,(i+j) mod p)).  Guests with pads are added on top.
A 'pad pair' is an intra-proto-cluster pair (v,i)-(v,k) that is taken out of the hard clique and made soft (owned by a pad member);
the proto-cluster stays glued as long as its remaining hard graph is connected."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *

class Builder:
    def __init__(self):
        self.names = []          # vertex names
        self.idx = {}
        self.E = set()
        self.members = []
        self.soft = set()
        self.hard_edges = set()
    def v(self, name):
        if name not in self.idx:
            self.idx[name] = len(self.names); self.names.append(name)
        return self.idx[name]
    def member(self, a, b, c):
        a, b, c = self.v(a), self.v(b), self.v(c)
        ps = pairs_of((a, b, c))
        assert not (ps & self.soft), ('pair reuse', self.names[a], self.names[b], self.names[c])
        self.soft |= ps
        self.members.append((a, b, c))
        self.E.add(key(a, b)); self.E.add(key(b, c))
    def hard(self, a, b):
        a, b = self.v(a), self.v(b)
        self.hard_edges.add(key(a, b))
    def finish(self):
        # hard edges must not be soft
        for p in self.hard_edges:
            assert p not in self.soft, ('hard edge is soft', [self.names[x] for x in p])
        E = set(self.E) | set(self.hard_edges)
        return len(self.names), E, list(self.members)

def blow_up(B, thin_members, thin_hard_pairs, p, vertex_names, pads_used=None):
    """thin_members: list of (a,b,c) over names; adds p^2 Latin members per thin member."""
    for (a, b, c) in thin_members:
        for i in range(p):
            for j in range(p):
                B.member((a, i), (b, j), (c, (i + j) % p))

def glue_clusters(B, names, p, pad_pairs):
    """hard clique on each proto-cluster except pad pairs (given as set of frozenset({(v,i),(v,k)}))."""
    for v in names:
        for i in range(p):
            for k in range(i + 1, p):
                pr = frozenset([(v, i), (v, k)])
                if pr in pad_pairs: continue
                B.hard((v, i), (v, k))
