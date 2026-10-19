"""Audit 7: Prop 5 of wave6_ml (family N*(F)): construction, properties F1-F6, and the ordering gap head(e_1) = tail(e_0)."""
import sys, itertools, subprocess, os, random
import numpy as np, networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in
import a7_lib as L, a7_ml as M

def build(edges):
    """edges: ordered list of (tail, head); edges[0] = (t, '*').  returns A, H, names, clique-of-vertex"""
    N1 = len(edges); idx = {}; cnt = 0
    for i in range(N1):
        idx[i] = dict(a=cnt, b=cnt + 1, c=cnt + 2, u=cnt + 3); cnt += 4
    z, y = cnt, cnt + 1; n = cnt + 2
    A = np.zeros((n, n), int)
    K = {}
    for i, (t, h) in enumerate(edges):
        K.setdefault(t, []).append(idx[i]['a'])
        K.setdefault(h, []).extend([idx[i]['b'], idx[i]['c'], idx[i]['u']])
    for v, vs in K.items():
        for p, q in itertools.combinations(vs, 2): A[p][q] = A[q][p] = 1
    def e(p, q): A[p][q] = A[q][p] = 1
    for i in range(N1):
        e(idx[i]['a'], idx[i]['b'])
        if i + 1 < N1: e(idx[i]['b'], idx[i + 1]['u'])
    e(idx[N1 - 1]['b'], z)
    e(y, idx[0]['u']); e(y, idx[0]['c'])
    H = []
    for i in range(N1): H.append((idx[i]['a'], idx[i]['b'], idx[i]['c']))
    for i in range(N1):
        u2 = idx[i + 1]['u'] if i + 1 < N1 else z
        H.append((idx[i]['u'], idx[i]['b'], u2))
    return A, H, idx, K, z, y

C_SRC = r'''
#include <stdio.h>
#include <stdlib.h>
int main(int argc,char**argv){ int n; scanf("%d",&n); int adj[32]={0}; for(int i=0;i<n;i++) for(int j=0;j<n;j++){int x; scanf("%d",&x); if(x) adj[i]|=1<<j;}
 long N=1L<<n; int *f=malloc(N*sizeof(int)); int *inner=malloc(N*sizeof(int)); inner[0]=0; f[0]=0;
 for(long m=1;m<N;m++){ int low=__builtin_ctzl(m); long r=m&(m-1); inner[m]=inner[r]+__builtin_popcount(adj[low]&r); int s=__builtin_popcountl(m); f[m]=s*(s-1)/2-2*inner[m]; }
 int *dp=malloc(N*sizeof(int)); dp[0]=0;
 for(long m=1;m<N;m++){ long low=m&-m; long rest=m^low; long sub=rest; int best=1<<30; while(1){ long c=sub|low; int v=dp[m^c]+f[c]; if(v<best) best=v; if(sub==0) break; sub=(sub-1)&rest; } dp[m]=best; }
 int E=0; for(int i=0;i<n;i++) E+=__builtin_popcount(adj[i]); E/=2; printf("%d\n",E+dp[N-1]); return 0; }
'''
def opt_c(A, tmpdir):
    src = os.path.join(tmpdir, "a7opt.c"); exe = os.path.join(tmpdir, "a7opt")
    if not os.path.exists(exe):
        open(src, "w").write(C_SRC); subprocess.check_call(["gcc", "-O2", "-o", exe, src])
    n = len(A); inp = str(n) + "\n" + "\n".join(" ".join(map(str, r)) for r in A) + "\n"
    return int(subprocess.run([exe], input=inp, capture_output=True, text=True).stdout)

def check(name, edges, tmpdir, do_opt=True, maxF=6):
    A, H, idx, K, z, y = build(edges); n = len(A); N1 = len(edges)
    print("==", name, "n", n, "|H|", len(H))
    ok_pack = L.is_packing(A, H); ok_pm = L.is_pair_maximal(A, H)
    I = M.comps_of(H); cmax = max(sum(1 for t in H if v in t) for v in range(n))
    Ig = nx.Graph(); Ig.add_nodes_from(range(len(H)))
    for i, j in itertools.combinations(range(len(H)), 2):
        if set(H[i]) & set(H[j]): Ig.add_edge(i, j)
    print(" packing", ok_pack, "pair-maximal", ok_pm, "components of I(H)", len(I), "tw(I(H))<=", treewidth_min_fill_in(Ig)[0], "c", cmax,
          "caterpillar/forest", nx.is_forest(Ig))
    S = L.soft_of(n, H); hard = A.astype(bool) & ~S
    D = 2 * cmax; l = 1; B = 2 * D + 3 * l
    print(" D", D, "max hard degree", int(hard.sum(1).max()), "B", B, "BIG empty:", int(hard.sum(1).max()) <= B)
    # F2: the partition
    blocks = list(K.values()) + [[z]]
    blocks = [list(b) for b in blocks]
    for b in blocks:
        if idx[0]['u'] in b and '*' is not None: pass
    # y into K_*^+ : the block that contains u_0
    for b in blocks:
        if idx[0]['u'] in b: b.append(y)
    lab = [-1] * n
    for k, b in enumerate(blocks):
        for v in b: lab[v] = k
    assert min(lab) >= 0
    cost = sum(1 for i, j in L.prs(n) if (lab[i] == lab[j]) != bool(A[i][j]))
    print(" cost(partition F2) =", cost, " |H|+1 =", len(H) + 1)
    if do_opt and n <= 26:
        print(" exact OPT =", opt_c(A, tmpdir))
    # F4 minor
    G2 = nx.Graph(); G2.add_nodes_from(range(n))
    for i, j in L.prs(n):
        if A[i][j] or S[i, j]: G2.add_edge(i, j)
    Q = nx.quotient_graph(G2, [set(b) for b in K.values()] + [{z}, {y}])
    print(" tw(Gamma') <= (min-fill)", treewidth_min_fill_in(G2)[0])
    # F5/F6
    cm = M.comps_of(H); assert len(cm) == 1
    W, sols = M.tight_solutions(A, H, cm[0])
    good = [1 for S_, B_ in sols if M.R_ok(A, B_, None, W)]
    print(" T tight:", len(sols) > 0, " T reducible:", len(good) > 0)
    # all pair-closed subfamilies F, |F| <= maxF
    pair_of = {}
    for i, t in enumerate(H):
        for p in L.p3p(t): pair_of[p] = i
    reducible = []; tested = 0
    for r in range(1, maxF + 1):
        for F in itertools.combinations(range(len(H)), r):
            Wf = sorted(set(v for i in F for v in H[i]))
            # pair-closed: no member outside F has a pair inside W_F
            if any(pair_of[p] not in F for p in itertools.combinations(Wf, 2) if p in pair_of): continue
            tested += 1
            W_, sols_ = M.tight_solutions(A, [H[i] for i in F] and H, list(F))
            if any(M.R_ok(A, B_, None, W_) for S_, B_ in sols_): reducible.append(F)
    print(" pair-closed subfamilies tested (|F|<=%d):" % maxF, tested, " reducible ones:", len(reducible), reducible[:5])
    return A, H

if __name__ == "__main__":
    tmp = sys.argv[1]
    # triangle, vertices 0,1,2 ; edges 0->1,1->2,2->0 ; e0 = (0,*)
    e0 = (0, '*')
    check("triangle, order e0,(0>1),(1>2),(2>0) [head(e1)=1 != tail(e0)]", [e0, (0, 1), (1, 2), (2, 0)], tmp)
    check("triangle, order e0,(2>0),(0>1),(1>2) [head(e1)=0 == tail(e0): GAP CASE]", [e0, (2, 0), (0, 1), (1, 2)], tmp)
