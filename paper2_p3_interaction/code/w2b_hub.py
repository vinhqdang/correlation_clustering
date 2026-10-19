"""Hub lemma check (deletion hubs / vertex-cover gadget are bounded).
Claim HUB: if X is a set of hard neighbours of u (ux in E, ux unpacked) that are pairwise non-adjacent
with every pair xx' unpacked, then the pairs of (x_1,u,x_2),(x_3,u,x_4),.. are all-hard, induced P3s of G,
modification-disjoint from each other and from H, so OPT >= |H| + floor(|X|/2).  Checked by brute force
for the inequality, using the maximum such X over all u.
Also prints the largest 'excess forced by a hub' seen (OPT-|H|) to show the bound is attained."""
import itertools, random
from check_csp_view import pairs, partitions
from w2b_sm import gen, rpack, opt
fs = frozenset
rng = random.Random(21); cnt = tight = 0
for _ in range(500):
    n = rng.randint(5, 8); E = gen(rng, n); H = rpack(rng, n, E, rng.random() < .5)
    packed = {q for P in H for q in pairs(P)}
    o = opt(n, E); ex = o - len(H)
    best = 0
    for u in range(n):
        Nh = [x for x in range(n) if fs((u, x)) in E and fs((u, x)) not in packed]
        for r in range(len(Nh), 1, -1):
            found = False
            for X in itertools.combinations(Nh, r):
                if all(fs((a, b)) not in E and fs((a, b)) not in packed for a, b in itertools.combinations(X, 2)):
                    found = True; break
            if found: best = max(best, r); break
    assert best // 2 <= ex, ("HUB", n, E, H)
    cnt += 1; tight += (best // 2 == ex and ex > 0)
print("HUB ok on", cnt, "instances; tight (floor(|X|/2)=excess>0) in", tight)
