"""Random check of Thm PL semantics: CLASS-OPT >= OPT always (soundness); CLASS-OPT = OPT whenever some optimum has hard edits inside U_loc (completeness on the class);
count how often the promise holds for random pair-maximal H with its min-fill decomposition."""
import random, sys, collections, time
from w8a_lib import *
from w8a_pl import class_opt, decomposition, uloc
rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2])
st = collections.Counter(); t0 = time.time()
for it in range(N):
    n = rng.randint(7, 10); E = [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < rng.choice([0.2, 0.3])]
    g = G8(n, E)
    if not g.p3s() or len(g.E) < 3: continue
    H = []; used = set(); P = g.p3s()[:]; rng.shuffle(P)
    for t in P:
        if not any(p in used for p in mp(t)): H.append(t); used.update(mp(t))
    opt, sols = g.opt(limit=500); l = opt - len(H)
    if l > 2: continue
    bags = decomposition(H); cand = set(uloc(g, H, bags))
    promise = any(set(hard_edits(H, edits_of_partition(g, part_from_mask(s, n)))) <= cand for s in sols)
    co, U, _ = class_opt(g, H, l, bags)
    sound = (co is None or co >= opt)
    complete = (co == opt)
    st[('l=%d' % l, 'promise=%s' % promise, 'sound=%s' % sound, 'CLASS-OPT=OPT: %s' % complete)] += 1
print('time %.0fs' % (time.time() - t0))
for k, v in sorted(st.items(), key=str): print('  ', k, v)
