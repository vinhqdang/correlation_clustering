"""T-B with slack: budget l = OPT-|H|+e.  Whole-component reduction of every tight reducible component (any size).
Checks: (1) OPT(G reduced) = OPT(G) - sum p_j over reduced comps; (2) #tight irreducible <= 2 l'; (3) #remaining comps <= l + l' <= 2 l;
(4) l=0 (e=0, OPT=|H|) => everything reducible."""
import sys, random
from w6a_lib import *
from w6a_rule import *
def main():
    seed, k, a, b = map(int, sys.argv[1:5]); rng = random.Random(seed)
    st = dict(inst=0, comps=0, reduced=0, opt_fail=0, k3_fail=0, r_fail=0, l0_fail=0, l0=0, maxrem=0, rem_gt_0=0)
    for it in range(k):
        n = rng.randint(a, b); adj = random_instance(n, rng); H = maximal_pair_packing(adj, rng)
        if not H: continue
        o = opt(adj); base = o - len(H)
        comps = I_components(H)
        info = []
        for T in comps:
            W = members_vertices(H, T); tau = opt(induced(adj, W))
            if tau > len(T): info.append((T, W, None, tau)); continue
            Ss = [S for S in tight_choices(adj, H, T, W) if R_ok(adj, S, W)]
            info.append((T, W, Ss[0] if Ss else False, tau))
        for e in (0, 1, 2):
            l = base + e; st["inst"] += 1
            lp = l - sum(tau - len(T) for T, W, S, tau in info)
            A = sum(1 for T, W, S, tau in info if S is False); B = sum(1 for T, W, S, tau in info if S is None)
            st["comps"] += len(info)
            if A > 2 * lp: st["k3_fail"] += 1
            if A + B > l + lp: st["r_fail"] += 1
            if base == 0 and e == 0:
                st["l0"] += 1
                if A or B: st["l0_fail"] += 1
            if e == 0:
                g2 = adj; red = 0
                for T, W, S, tau in info:
                    if S: g2 = flip(g2, S); red += len(T); st["reduced"] += 1
                if opt(g2) != o - red: st["opt_fail"] += 1; print("OPT FAIL", sorted(edges_of(adj)), H)
                st["maxrem"] = max(st["maxrem"], A + B); st["rem_gt_0"] += (A + B > 0)
    print(st)
main()
