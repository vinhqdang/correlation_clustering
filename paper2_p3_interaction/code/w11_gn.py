"""w11_gn.py -- the persistent-bit chain G_N (definition from notes/wave9_corebound.md Sec 5.1) and its certificates, new code."""
import sys, time
from w11_lib import *

def build_GN(N, pendants=True):
    vid = {}
    def V(j, name):
        if name == 'a' and j > 0:
            return V(j - 1, 'b')
        if (j, name) not in vid:
            vid[(j, name)] = len(vid)
        return vid[(j, name)]
    H = []; E = set()
    for j in range(N):
        p, q, r, c, b, a = (V(j, x) for x in 'pqrcba')
        H += [(q, p, a), (q, r, b), (p, c, r)]
        E |= {key(q, p), key(p, a), key(q, r), key(r, b), key(p, c), key(c, r), key(c, q)}
    nv = len(vid)
    T = []
    if pendants:
        for j in range(N):
            w = nv; nv += 1
            E.add(key(w, V(j, 'p'))); T += [key(w, V(j, 'p')), key(w, V(j, 'r'))]
        w = nv; nv += 1
        E.add(key(w, V(0, 'a'))); T += [key(w, V(0, 'a')), key(w, V(0, 'p'))]
        w = nv; nv += 1
        E.add(key(w, V(N - 1, 'b'))); T += [key(w, V(N - 1, 'b')), key(w, V(N - 1, 'r'))]
    return Inst(nv, E, H), T, V

if __name__ == '__main__':
    for N in range(1, int(sys.argv[1]) + 1):
        I, T, V = build_GN(N)
        t0 = time.time()
        O = FPrime(I, len(I.H))
        Tn = [e for e in T]
        print('N', N, 'n', I.n, 'valid', I.valid(), 'pm', I.pair_maximal(), 'c', I.cmax(), 'softtri', len(I.soft_triangles()),
              'F(T_N) feasible at |H|?', O.feasible(Tn) is not None,
              'all F(T-e) feasible?', all(O.feasible([x for x in Tn if x != e]) is not None for e in Tn),
              'remote ⊇ T', set(Tn) <= set(I.remote()))
        sz, Tm = min_certificate(I, 0)
        print('   min certificate l=0:', sz, [ (x,y) for x,y in Tm], round(time.time() - t0, 1), 's', flush=True)
        if N <= 3:
            print('   OPT - |H| =', opt(I) - len(I.H) if I.n <= 16 else 'n/a')
