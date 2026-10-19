"""exhaustive chain search: path of m members with a given sharing pattern, all graphs with <= kx extra edges (non-member pairs).
Keep: pair-maximal, OPT = |H| (tight solution exists), T reducible as a whole, NO reducible proper connected subfamily."""
import sys, itertools, collections
from w6a_lib import *
from w6a_rule import *
def pattern_members(m, pat):
    # pat[i-1] = (role_prev, role_new): role in {'e','c'} ; member tuple (x,y,z): y centre
    mem = [(0, 1, 2)]; nv = 3; used = [set()]
    for i in range(1, m):
        rp, rn = pat[i - 1]
        prev = mem[-1]
        pos = {'e': [0, 2], 'c': [1]}[rp]
        pos = [p for p in pos if prev[p] not in used[-1]]
        if not pos: return None
        v = prev[pos[0]]; used[-1].add(v)
        a, b = nv, nv + 1; nv += 2
        t = [a, b]
        t.insert(0 if rn == 'e' else 1, v)   # 'e': v is endpoint x ; 'c': v is centre
        mem.append(tuple(t)); used.append({v})
    return mem, nv
def run(m, kx, pat):
    r = pattern_members(m, pat)
    if r is None: return []
    mem, n = r
    forced = set(); E0 = set()
    for x, y, z in mem:
        E0 |= {key(x, y), key(y, z)}; forced |= {key(x, y), key(y, z), key(x, z)}
    free = [key(a, b) for a in range(n) for b in range(a + 1, n) if key(a, b) not in forced]
    out = []
    T = list(range(m)); W = members_vertices(mem, T)
    subs = [F for F in conn_subfamilies(mem, m - 1)]
    for k in range(0, kx + 1):
        for ex in itertools.combinations(free, k):
            adj = from_edges(n, E0 | set(ex))
            if not is_pair_maximal(adj, mem): continue
            ch = tight_choices(adj, mem, T, W)
            if not ch: continue
            if not any(R_ok(adj, S, W) for S in ch): continue
            bad = False
            for F in subs:
                W2 = members_vertices(mem, F)
                if pair_closed(mem, F, W2) and any(R_ok(adj, S, W2) for S in tight_choices(adj, mem, F, W2)): bad = True; break
            if bad: continue
            out.append((n, sorted(E0 | set(ex)), mem))
    return out
if __name__ == "__main__":
    m = int(sys.argv[1]); kx = int(sys.argv[2])
    pats = list(itertools.product(itertools.product('ec', 'ec'), repeat=m - 1))
    seen = set()
    for pat in pats:
        res = run(m, kx, pat)
        print(m, kx, pat, len(res), flush=True)
        for r in res[:2]: print("   ", r, flush=True)
