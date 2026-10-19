"""Wave 6: tests of the subfamily / component-level Rule 6.2 (BFK) on pair-maximal packings.
 T-A  soundness: tight S_F with R-condition  =>  OPT(G^S) = OPT(G) - |F|   (every such S, every connected pair-closed F, |F|<=4)
 T-B  component level: counts (K3-type) #tight-irreducible <= 2 l',  #components after reduction <= l + l'
 T-C  how often single-member rules suffice; examples of components reducible only as a whole
Usage: python3 w6a_rule.py seed ninst nmin nmax
"""
import sys, random, itertools
from w6a_lib import *

def members_vertices(H, F):
    W = 0
    for i in F:
        for v in H[i]: W |= 1 << v
    return W

def conn_subfamilies(H, maxsize):
    m = len(H); Ie = Iedges(H)
    nb = [set() for _ in range(m)]
    for i, j in Ie: nb[i].add(j); nb[j].add(i)
    seen = set(); out = []
    def grow(cur, ext, forb):
        out.append(tuple(sorted(cur)))
        if len(cur) == maxsize: return
        for x in sorted(ext):
            if x in forb: continue
            forb2 = forb | {x}
            # standard enumeration of connected subsets
        return
    # simple BFS enumeration
    level = {frozenset([i]) for i in range(m)}
    allf = set(level)
    for _ in range(maxsize - 1):
        nxt = set()
        for F in level:
            ext = set()
            for i in F: ext |= nb[i]
            ext -= F
            for x in ext: nxt.add(F | {x})
        allf |= nxt; level = nxt
    return [tuple(sorted(F)) for F in allf]

def pair_closed(H, F, W):
    Fs = set(F)
    for j, t in enumerate(H):
        if j in Fs: continue
        if sum(1 for v in t if W >> v & 1) >= 2: return False
    return True

def tight_choices(adj, H, F, W):
    """all S (one pair edit per member of F) such that G[W] xor S is a cluster graph"""
    res = []
    opts = [list(pairs3(H[i])) for i in F]
    for ch in itertools.product(*opts):
        if len(set(ch)) < len(ch): continue
        a2 = flip(adj, ch)
        ok = True
        for v in range(len(adj)):
            if W >> v & 1:
                nv = a2[v] & W
                x = nv
                while x:
                    u = (x & -x).bit_length() - 1; x &= x - 1
                    if ((a2[u] & W) | 1 << u) != (nv | 1 << v): ok = False; break
                if not ok: break
        if ok: res.append(ch)
    return res

def R_ok(adj, S, W):
    a2 = flip(adj, S); n = len(adj)
    clusters = set()
    for v in range(n):
        if W >> v & 1: clusters.add((a2[v] & W) | 1 << v)
    for y in range(n):
        if W >> y & 1: continue
        A = adj[y] & W
        if A and A not in clusters: return False
    return True

def induced(adj, W):
    vs = [v for v in range(len(adj)) if W >> v & 1]
    idx = {v: i for i, v in enumerate(vs)}
    E = {(idx[a], idx[b]) for a in vs for b in vs if a < b and adj[a] >> b & 1}
    return from_edges(len(vs), E)

def check_soundness(adj, H, base_opt, maxsize, stats, rng):
    for F in conn_subfamilies(H, maxsize):
        W = members_vertices(H, F)
        if not pair_closed(H, F, W): continue
        for S in tight_choices(adj, H, F, W):
            if R_ok(adj, S, W):
                stats["applic"] += 1
                o2 = opt(flip(adj, S))
                stats["sound_fail"] += (o2 != base_opt - len(F))
                if o2 != base_opt - len(F): print("SOUNDNESS FAIL", adj, H, F, S)

def run(seed, ninst, nmin, nmax):
    rng = random.Random(seed)
    st = dict(inst=0, applic=0, sound_fail=0, comp_tot=0, comp_tight=0, comp_nontight=0, comp_red=0,
              k3_fail=0, r_fail=0, multi_comp=0, whole_only=0, whole_only_ex=None, nonpm=0)
    for it in range(ninst):
        n = rng.randint(nmin, nmax)
        adj = random_instance(n, rng)
        H = maximal_pair_packing(adj, rng)
        if not H: continue
        assert is_pair_maximal(adj, H)
        o = opt(adj); l = o - len(H)
        if l < 0: print("BOUND FAIL"); continue
        st["inst"] += 1
        check_soundness(adj, H, o, 3, st, rng)
        # component level
        comps = I_components(H)
        irr = 0; nt = 0; taus = 0; tight_irr = 0
        for T in comps:
            W = members_vertices(H, T)
            tau = opt(induced(adj, W))
            taus += tau - len(T)
            if tau > len(T): nt += 1; continue
            red = any(R_ok(adj, S, W) for S in tight_choices(adj, H, T, W)) if pair_closed(H, T, W) else False
            # components are always pair closed (distinct components share no vertex)
            if not red: tight_irr += 1
            else:
                st["comp_red"] += 1
        lp = l - taus   # excess l' = l - sum (tau_j - p_j)  (l = OPT-|H|)
        st["comp_tot"] += len(comps); st["comp_nontight"] += nt
        if tight_irr > 2 * lp: st["k3_fail"] += 1; print("K3 FAIL", adj, H)
        if tight_irr + nt > l + lp: st["r_fail"] += 1
        # whole-component-only reducibility: component reducible, but no connected subfamily of size<=2 reducible
        for T in comps:
            if len(T) < 3: continue
            W = members_vertices(H, T)
            if opt(induced(adj, W)) > len(T): continue
            if not any(R_ok(adj, S, W) for S in tight_choices(adj, H, T, W)): continue
            sub_red = False
            for F in conn_subfamilies([H[i] for i in T], 2):
                FF = tuple(T[i] for i in F)
                W2 = members_vertices(H, FF)
                if pair_closed(H, FF, W2) and any(R_ok(adj, S, W2) for S in tight_choices(adj, H, FF, W2)): sub_red = True; break
            if not sub_red:
                st["whole_only"] += 1
                if st["whole_only_ex"] is None: st["whole_only_ex"] = (n, sorted(edges_of(adj)), H, T, l)
    print(st)

if __name__ == "__main__":
    s, k, a, b = map(int, sys.argv[1:5]); run(s, k, a, b)
