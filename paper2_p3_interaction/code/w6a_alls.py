"""N*_path(F): test ALL (not nec. connected) subfamilies F, |F|<=smax, for Rule-6.2 reducibility (tight S_F, R w.r.t. G, pair-closed)."""
import sys, itertools
from w6a_lib import *
from w6a_rule import *
from w6a_csp import tight_choices_dfs
from w6a_trig3 import build_path
def main():
    tri = [(0, 1), (1, 2), (2, 0)]
    K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for name, ed, smax in [("triangle+pendant", tri, 8), ("K4+pendant", K4, 5)]:
        n, adj, H, cl, y, idx, ed2 = build_path(ed)
        tot = 0; red = []; pc_ = 0
        for k in range(1, smax + 1):
            for F in itertools.combinations(range(len(H)), k):
                W = members_vertices(H, F)
                if not pair_closed(H, F, W): continue
                pc_ += 1
                if any(R_ok(adj, S, W) for S in tight_choices_dfs(adj, H, F, W)): red.append(F)
        print(name, "members", len(H), "pair-closed subfamilies tested", pc_, "reducible:", len(red), red[:5], flush=True)
main()
