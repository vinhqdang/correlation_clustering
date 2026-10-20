"""Independent verification of a synthesized chain instance (JSON from w8c_per.py).
Checks: valid packing; anchors clean; ZE exists; for every subset S of anchors: ZE with S singletons exists?  -> face table;
exact OPT by ILP (HiGHS) of the graph with stubs on S (a stub = new pendant vertex adjacent to the anchor only)."""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *

def main(path, ilp=True):
    d = json.load(open(path))
    n = d['n']; E = set(tuple(e) for e in d['E']); members = [tuple(m) for m in d['members']]; A = d['anchors']
    print('n', n, '|H|', len(members), 'anchors', A)
    print('valid', check_instance(n, E, members))
    soft = set()
    for m in members: soft |= pairs_of(m)
    for x in A:
        hard = [v for v in range(n) if v != x and key(x, v) in E and key(x, v) not in soft]
        cnt = sum(1 for m in members if x in m)
        print(' anchor', x, d['vname'][x], 'hard edges', hard, 'members through', cnt, 'role', [ 'end' if m[1]!=x else 'centre' for m in members if x in m])
    # maximum c
    c = max(sum(1 for m in members if v in m) for v in range(n)); print('c =', c)
    # tw of I(H)
    import itertools as it
    from w6c_lib import tw_heur, ih_graph
    print('tw(I(H)) <= (min-fill)', tw_heur(ih_graph(members)))
    faces = {}
    for r in range(len(A) + 1):
        for S in it.combinations(A, r):
            faces[S] = ze_solve(n, E, members, singletons=list(S)) is not None
    for S, f in faces.items():
        if len(S) >= len(A) - 1: print(' lock', [d['vname'][x][0] for x in S], 'ZE feasible:', f)
    minnon = [S for S, f in faces.items() if not f and all(faces[T] for T in it.combinations(S, len(S) - 1))]
    print('minimal non-faces (cores):', [[d['vname'][x][0] for x in S] for S in minnon])
    if ilp:
        for S in [tuple(A)] + [tuple(x for x in A if x != a) for a in A]:
            E2 = set(E); n2 = n
            for x in S:
                E2.add(key(x, n2)); n2 += 1
            opt = ilp_opt_ok(n2, E2, members)
            print(' ILP OPT with stubs on', [d['vname'][x][0] for x in S], '=', opt, ' |H| =', len(members))

if __name__ == '__main__':
    main(sys.argv[1], ilp=(len(sys.argv) < 3 or sys.argv[2] != 'noilp'))
