# Wave 10 (prover C'): W[1]-hardness of CEaMP in l + tw(I(H))  --  comparator found, reduction assembled, proof written

Labels: PROVED = proof in this file (a finite computation that is part of a proof is named and marked "by computer");
VERIFIED-BY-COMPUTER = exact finite check of the stated instances only (script named); HEURISTIC; OPEN.
Code: `code/w10_*.py` (library `w10_lib.py`, reduction `w10_red.py`; uses `w8c_lib.py`, `w8c_oh.py`, `w8c_thick.py`, `w8c_chainbuild.py` of wave 8).
No git, `paper/` untouched. Notation as in wave6/8/9: soft pair = pair inside a member, hard pair = other pair, member (a,b,c) has centre b,
states L/R/M = delete ab / delete bc / insert ac, ZE = zero excess (every member edited exactly once, nothing else), ex(S) = |S| - |H|.

## 0. Verdict

1. **Main result (Theorem H).** CEaMP is W[1]-hard parameterised by `l + tw(I(H))`, and even by `l + c + pw(I(H))`; the hard instances are pair-maximal.
   Reduction from Grid Tiling with <= (Marx 2007; Cygan et al., *Parameterized Algorithms*, Thm 14.28): instance (k,n,S) -> (G,H,l) with
   `l = k^2`, `c = 4p+2`, `pw(I(H)) <= 22 k^2 p^2` (measured per chain 651 = 10.2 p^2 at p = 8, 1014 at p = 10, 1457 at p = 12; totals 1303 for k=1, 5199 for k=2, 11719 for k=3 at p = 8; independent of n), where `p` is an even number `>= k^2 + 5`; size poly(k,n).
   So all three parameters are functions of k only. (Parameter values are polynomial in k: `l=k^2`, `c=O(k^2)`, `pw=O(k^6)`.)
   **Status of the proof:** every step is PROVED; the proofs use three finite computations (chain certificate, profile table, explicit-cluster check),
   each reproduced by a script and valid for all m, and the whole construction is VERIFIED-BY-COMPUTER end-to-end on tiny instances (Sec. 6).
   The only element that is not a full written proof with a symbolic argument is the *finite-state certificate for the chain* (Lemma 1(b)), which is a
   transfer-matrix computation (explicit 6x6 integer matrix, char. poly `x^4 (x-1)^2`); I label it PROVED-by-computer.
2. **The missing comparator exists and is simple**: not a gadget between two chains of the old kind, but an *OR-switch guest*
   (Sec. 2.4): one new vertex `y` and one member `(z_A, y, z_B)` where `z_A,z_B` are "permanent guests" hosted in a chain-A / chain-B cluster.
   `y` must sit with exactly one of `z_A, z_B`; sitting with `z_A` is possible iff the chain-A position lies in a threshold set X, with `z_B` iff the chain-B position lies in a threshold set Y.
   It realises the clause `(a in X) or (b in Y)`; the comparator `a <= b` at block granularity is the conjunction of `n-1` such clauses with `X = {a <= t}`, `Y = {b > t}`.
   VERIFIED (thin and thick): the 10x10 table of position pairs (a,b) is feasible exactly on `(a<=3 or b>=6)` minus the permanently forbidden positions (Sec. 6.1).
3. **Soundness is a counting argument, not a case analysis** (Sec. 5): `ex >= max(q, 2k^2 - q) >= k^2` where `q` = number of edited cross-anchor pairs; at `ex = k^2` the solution is forced
   to be "one cut tunnel per cell + ZE otherwise". This is the "excess >= #tracks/2" accounting of wave 8, now with a complete proof including the glue (needs `p >= l+5`) and
   the comparator guests (which lie in no chain region).
4. Consistency with wave 9: Theorem A / Helly-2 forbid hardness on soft-triangle-free instances; the instances here are full of soft triangles (Latin members of the blow-up)
   and have large c. **Whether hardness survives with c bounded (e.g. c = 34 constant, or c = 2) is OPEN**; in this proof p (hence c) must grow with l because the glue must survive `l` edits (Sec. 5, Lemma 5).

## 1. Conventions and tools

* `OPT(G)`, solutions, decomposition lemma (paper prelim): every solution S has `x_i = S cap pairs(P_i) != {}` and `|S| = sum |x_i| + |U|`, U = hard edits; so `ex(S) >= |U|`.
* ZE (zero excess) partition: exactly one edit per member, no hard edit; per member state L/R/M.
* Thin chain `OH_m` (wave 8, `w8c_oh.py`, m >= 7): vertices `s_0..s_{2m}`, `e_0..e_m`; members `P_j=(s_{2j},s_{2j+1},s_{2j+2})` (spine) and interior members
  `G_j=(e_j,s_{2j+1},e_{j+2})`, `F_j=(e_j,s_{2j+2},e_{j+1})` plus fixed left/right caps; four hard edges; every other pair is a hard non-edge. *Position* a in {0..m-1}: the spine word `R^a M L^{m-1-a}`.
* Thick: each thin vertex v -> *super-vertex* `K_v = {(v,i): i<p}`, p even; each thin member (a,b,c) -> p^2 Latin members `((a,i),(b,j),(c,i+j mod p))`;
  each thin hard pair -> all p^2 thick pairs; inside `K_v` all pairs are hard edges except pairs taken as *pad pairs* (below).
  Every pair between the super-vertices of one thin member is covered by exactly one Latin member (soft).
* **Guest** = a new single vertex y with a set `W_y` of super-vertices it may share a cluster with. For `v in W_y`: take a perfect matching `M` of `K_p` (one unused 1-factor of the round-robin 1-factorisation per (guest, v)),
  add the p/2 *pad members* `(y,(v,i),(v,k))` with centre `(v,i)`, `(i,k) in M` (the pair `(v,i)(v,k)` is a soft edge). All other pairs of y are hard non-edges.
  Super-vertices are glued by the hard clique minus all pad pairs (connected: checked by the builder for every v).
* *Permanent guest* z at host h: members `((h,i),z,(h,k))` for a matching (centre z; the pair `(h,i)(h,k)` is a soft non-edge) and pads to `W_z`.
* Tools: `ze_solve_fast`, `ze_atm` (SAT, CaDiCaL), `exact_decide` (exact OPT decision by lazy-triple SAT, Sec. 6.5), HiGHS not needed.

## 2. Construction

### 2.1 Chains, positions, groups
All chains are copies of the thick `OH_m` with `m = 2U`, `U = 1 + n(n+1)` (for n >= 2; U >= 7). *Group* `u in {0..U-1}` = positions `{2u, 2u+1}`; its *host* is `h_u = s_{4u+2}`.
Block `alpha in [n]` has the slot groups `u = (alpha-1)(n+1)+s`, `s in [n]`, followed by the separator group `u = alpha(n+1)`; group 0 is an end pad (the last separator is the end pad on the right).

### 2.2 Explicit clusters and classes (Lemma 1/2 below)
For host `h_u` the positions fall into four classes with the same host cluster: `A0 = {a <= 2u-1}`, `A1 = {2u}`, `A2 = {2u+1}`, `A3 = {a >= 2u+2}`.
Their partner sets `P(A_i)` (cluster of `h_u` minus `h_u`): `P(A0) = {s_{4u+1}, eps(2u)}`, `P(A1) = {s_{4u},s_{4u+1},eps(2u)}`, `P(A2)={s_{4u+3},s_{4u+4},eps(2u+1)}`, `P(A3)={s_{4u+3},eps(2u+1)}`
(eps = the e-label of the cell, with the left/right cap exceptions of `explicit_partition` in `w10_chain.py`). Facts used: `P(A0) subset P(A1)`, `P(A3) subset P(A2)`, `s_{4u} in P(A1) \ P(A0)`, `s_{4u+4} in P(A2) \ P(A3)`, `P(A0) cap P(A2) = P(A0) cap P(A3) = empty`.
Four allowed sets and their *effective forbidden sets* (positions whose host cluster is not inside `W_y`):

| name | W | forbidden positions |
|---|---|---|
| pad | host + P(A0) + P(A3) | `{2u,2u+1}` |
| X | host + P(A0) | `{a >= 2u}` (so the guest may join iff `a <= 2u-1`) |
| Y | host + P(A3) | `{a <= 2u+1}` (iff `a >= 2u+2`) |
| all | host + everything | none |

### 2.3 Guests per chain
For every group u of chain c: a permanent guest `z_{c,u}` at host `h_u`. If u is a *tile group* (a slot group of a tile of the cell, see 2.5): `W_z = all` and, in addition, a *lock guest* `y_{c,u}` with `W_y = pad`,
an *anchor* `x_{c,u}`, and the member `(x_{c,u}, y_{c,u}, z_{c,u})`. Otherwise (slot of a non-tile, separator, end pad): `W_z = pad` (a *permanent lock*: z cannot be in the host cluster at the two positions of its group).

### 2.4 The comparator (OR-switch guest)
For two chains `c1, c2` and a separator `alpha < n` (group `u_alpha`, host `h`): a guest `y` with `W_y = X(c1) + Y(c2)` (pads to the super-vertices of `X` in chain `c1` and of `Y` in chain `c2`) and the member `(z_{c1,u_alpha}, y, z_{c2,u_alpha})`.
Semantics (Prop. 4): y is with `z_{c1}` (needs `a_{c1} <= 2u_alpha - 1`) or with `z_{c2}` (needs `a_{c2} >= 2u_alpha + 2`), never both, never none.
Note that `z_{c,u_alpha}` serves up to two comparators (left neighbour as second chain, right neighbour as first chain); their position sets are disjoint, so the two guests never meet in a cluster.

### 2.5 Grid Tiling instance
`GT(<=)`: k, n, `S_ij subset [n]^2` (nonempty), find `s_ij=(a_ij,b_ij) in S_ij` with `a_ij <= a_{i,j+1}` and `b_ij <= b_{i+1,j}`.
Chains `A_ij, B_ij` for every cell. Tile `(a,b) in S_ij` has its anchor at slot group `(alpha=a, s=b)` of `A_ij` and at slot group `(alpha=b, s=a)` of `B_ij` (the two orders (a,b) and (b,a)); all other slots are permanent locks.
*Tunnel*: hard edge between the two anchors of the same tile. Comparators: `(A_ij, A_{i,j+1}, alpha)` and `(B_ij, B_{i+1,j}, alpha)` for `alpha = 1..n-1`. Budget `l = k^2`.
`I_tile` below denotes this instance. Anchors have no other hard edge; every tunnel joins `A_ij` and `B_ij` of the same cell.

## 3. Lemmas about the chain and the guests

**Lemma 1 (chain; PROVED, finite certificate by computer).** For every m >= 7 the ZE partitions of `OH_m` are exactly `P_a`, a = 0..m-1, given by the explicit formula `explicit_partition(m,a)` of `w10_chain.py`
(clusters: cell j of the spine contributes `{s_{2j},s_{2j+1}}` if j < a, `{s_{2j},s_{2j+1},s_{2j+2}}` if j = a, `{s_{2j+1},s_{2j+2}}` if j > a, together with the e-vertex of the cell; two cap exceptions).
(a) Each `P_a` is ZE: all members have one edit, hard pairs respected. This is a conjunction of constraints each involving vertices of at most 3 consecutive cells; windows in the interior are isomorphic,
left-cap windows do not depend on m, right-cap windows depend only on the distance to the right end; hence the check for m = 7..60 (all a; `w10_chain.py`: 0 invalid) covers all m.
(b) There are no other ZE partitions: the cell-by-cell dynamic programme (state = block structure of the vertices of the last two cells, `w8c_transfer2.step_states`) has exactly 6 states at every interior cell and
an interior transfer matrix M that is **identical for all cells 7..19 of m = 26** (and, by translation invariance of the member list, for every interior cell of every m):
`M = [[0,0,0,0,1,1],[0,0,0,0,0,0],[1,1,1,1,0,0],[0,0,0,0,0,0],[0,0,0,0,1,1],[0,0,0,0,0,0]]`, characteristic polynomial `x^4 (x-1)^2`, Jordan block of size 2 at eigenvalue 1.
Hence the number of ZE partitions `N_m = u^T M^{m-c0} v` is, for m >= m1, an affine function of m; `N_m = m` for m = 7..40 (`w8c_transfer.py`, `w8c_transfer2.py`), so `N_m = m` for all m >= 7.
Together with (a) the m partitions `P_a` are all of them. (Run `python3 w10_chain.py` and `python3 w10_chain.py transfer`.)

**Lemma 2 (profile / guest lemma; PROVED from Lemma 1(a)).** For the host `h_u` of every group the classes `A0..A3` of Sec. 2.2 have the stated partner sets; the allowed sets pad/X/Y/all have exactly the effective forbidden sets of the table in Sec. 2.2.
Proof: read the host cluster off the formula of Lemma 1(a): `s_{4u+2}` is the `c`-vertex of cell `2u` and the `a`-vertex of cell `2u+1`; cell `2u` is in state L for a <= 2u-1 (cluster `{s_{4u+1},s_{4u+2}}+e`),
M for a = 2u (`{s_{4u},s_{4u+1},s_{4u+2}}+e`), and cell 2u+1 is in state M for a = 2u+1, R for a >= 2u+2. The e-labels of cells 2u and 2u+1 are different (also at the caps). The membership claims are then immediate.
Checked by computer for all even m = 8..60 and all u (1,836 (m,u,type) cases, 0 failures: `w10_prof.py`), including the cap groups u = 0, U-1.

**Lemma 3 (blow-up; PROVED).** Let G be a thick instance (Sec. 1) and `G_thin` the *thin-guest model*: the thin chain(s) plus guests as single vertices, with the pair `(y,K_v)` for `v in W_y` *free* (no cost, no constraint other than transitivity),
the pair (z,host) a hard edge, all other pairs as in G. Then the ZE partitions of G are exactly the lifts of the ZE partitions of `G_thin`.
Proof. In a ZE partition no hard pair is edited, so the connected hard glue puts each `K_v` into one cluster. A thin hard pair gives p^2 hard pairs (edge: together, non-edge: apart).
The p^2 Latin members of a thin member see the same togetherness pattern of three super-vertices, so they are in the same state L/R/M and this is a thin state; conversely lifting a thin ZE partition gives one edit per Latin member.
A pad member `(y,(v,i),(v,k))` with `(v,i)~(v,k)` has exactly one edit in both cases: y together with `K_v` (state M: insert `y(v,k)`) or apart (state L: delete `y(v,i)`); it is never in state R. If `v notin W_y` the pair `y(v,i)` is a hard non-edge, so y is apart from `K_v`.
The host members `((h,i),z,(h,k))` force state M (since `(h,i)~(h,k)`), i.e. `z ~ K_h`. QED. (VERIFIED thick vs thin: Sec. 6.1, 6.2.)

**Proposition 4 (tunnel-release semantics TR; PROVED from Lemmas 1-3).** Let `C` be a set of tunnels containing exactly one tunnel of each chain (a tunnel serves its two chains), and `G' = I_tile - C` (the tunnel edges of C removed as pairs of G, others kept).
Then `G'` has a ZE partition iff the cut tiles `t_ij` (the tile of the cut tunnel in cell (i,j)) satisfy the Grid Tiling constraints.
Moreover in every ZE partition chain c is in a position `a_c` of the group of its cut tunnel.
Proof. (=>) Chain c: by Lemma 1 the restriction to c is some `P_{a_c}`. An uncut tunnel `x x'` is an unedited hard edge, so `x ~ x'`; any further member of the cluster of x would have to be soft to x and to x' (all other pairs
at x, x' are hard non-edges), but `Q(x) cap Q(x')` is empty, so the cluster is `{x,x'}` (Lemma 2 of wave 6). In the member `(x,y,z)` the state is therefore L, i.e. `y ~ z ~ K_host` (Lemma 3), so `K_host`'s cluster at position `a_c` is inside `W_y`:
`a_c notin {2u,2u+1}`. Permanent locks give the same for the other groups. Since the groups partition [0,m), `a_c` lies in the group of the cut tunnel.
For a separator `u_alpha` with comparator guest y between chains c1,c2: y is with `z_{c1}` or with `z_{c2}` (state L or R of `(z_{c1},y,z_{c2})`; M would put `z_{c1}` and `z_{c2}` into one cluster, impossible since they lie in clusters of different chains),
so `a_{c1} <= 2u_alpha - 1` or `a_{c2} >= 2u_alpha+2` (Lemma 2, table). Writing `beta(c)` for the block of the cut tile in c: `a_c <= 2u_alpha-1` iff `beta(c) <= alpha`, `a_c >= 2u_alpha+2` iff `beta(c) >= alpha+1`.
So for every alpha: `beta(c1) <= alpha or beta(c2) > alpha`, i.e. `beta(c1) <= beta(c2)`. For rows `beta(A_ij) = a_ij`, for columns `beta(B_ij) = b_ij` (the B-chain order is by (b,a)), which are the constraints of GT(<=).
(<=) Given a tiling, cut the tunnel of the tile `s_ij` in every cell; put chain c at position `2u` (u the group of the cut tunnel). Place guests: `z`, locked `y` (with `x x'` as clusters `{x,x'}`) join the host clusters (allowed by the table: the host cluster at `a_c` is inside `W` because `a_c` is outside the group of every locked anchor and every permanent lock);
the released anchor x (no hard edge) forms the cluster `{x,y}` (state R, y apart from z); each comparator guest joins `z_{c1}` if `beta(c1) <= alpha` (allowed by the table) and `z_{c2}` otherwise (then `beta(c2) >= alpha+1`).
Two guests never share a cluster: guests join a base cluster only through a permanent guest z, hosts `s_{4u+2}` and `s_{4u'+2}` are never in one cluster (positions differ by >= 4), and the two comparator guests on the same z have disjoint position sets.
All remaining pairs are as required. QED. (VERIFIED: Sec. 6.1-6.3.)

## 4. Parameters (Theorem P)

* `l = k^2`; `c = 4p+2` (every vertex `(v,i)` lies in p Latin members of each of at most 4 thin members plus 2 guest members; measured 34, 42, 50 for p = 8, 10, 12).
* `|V| = O(k^2 n^2 p)`, `|H| = O(k^2 n^2 p^2)`: polynomial.
* `tw(I(H)) <= pw(I(H)) <= 22 k^2 p^2` (measured, fits `10.2 p^2` per chain for p = 8,10,12 and any k, n): all members have a coordinate interval of length <= 2 cells (a member of cell j touches cells j..j+2; guests and comparators sit at the cell `2u+1` of their group, a comparator member joins two chains
at the *same* cell), tunnels are hard pairs and do not appear in `I(H)`, so sweeping the cells from left to right with all 2k^2 chains in parallel gives bags of size `<= 2k^2 * (members through one cell)`;
the maximum cover of the interval supergraph is 1303 (k=1) / 5199 (k=2) / 5207 (k=2, n=3) / 5199 (k=2, n=4) for p = 8 (`w10_tw.py`): independent of n, linear in the number of chains, proportional to p^2.
* `p`: even, `p >= l + u_max + 2` with `u_max = 3` (Lemma 5), i.e. `p = k^2 + 5` rounded up to even.

## 5. Soundness (Theorem S) and completeness (Theorem C)

**Lemma 5 (no super-vertex is split; PROVED).** Let `u_max` be the maximum number of pad/host matchings used at a super-vertex (measured 2 for k <= 2, 3 for k >= 3; by construction <= 3).
If `p >= l + u_max + 2` and S is a solution with `ex(S) <= l`, then every `K_v` lies in one cluster of `G xor S`; in particular no pair inside a super-vertex is edited.
Proof. Suppose `K_v` meets two clusters; let `X` be the part in one of them, `1 <= |X| <= p/2`, `Y` the rest. All pairs `X x Y` inside `K_v` that are edges must be deleted. Pairs of `K_v` that are soft are the pad/host pairs: at most `u_max` per vertex, and the host pairs are non-edges.
So at least `|X|(|Y| - u_max) >= p - 1 - u_max > l` pairs are *hard* edges that are edited, contradicting `|U| <= ex(S) <= l` (decomposition lemma). QED.

**Regions.** For a chain c let `R_c` = all its thick vertices and its guests `z,y,x` (not the comparator guests). The regions are pairwise disjoint. `T_c` = the tunnel partners of the anchors of c (they lie in the partner chain `c'`), `H_c` = the members inside `R_c`
(chain, permanent-guest, lock and pad members of its own guests). `H_cmp` = members containing a comparator guest (the pad members and the switch members), `H = union H_c + H_cmp`.
`Pi_c` = pairs inside `R_c cup T_c` except pairs with both ends in `T_c`. A pair lies in two `Pi`'s iff it joins an anchor of c to an anchor of the partner chain (a *cross-anchor pair*, `q` := number of edited ones); otherwise in at most one.

**Lemma 6 (region bound; PROVED from Lemmas 1-3).** Any solution S has `|S cap Pi_c| >= |H_c| + 1`. Proof: the members of `H_c` are pair-disjoint and each has an edit, so `|S cap Pi_c| >= |H_c|`; equality means one edit per member of `H_c` and *no other edit in Pi_c*.
Then the partition restricted to `R_c cup T_c` is ZE for the graph `G_c = G[R_c cup T_c]` (stubs: each partner anchor of `T_c` has as only pairs the tunnel edge to its anchor in c, all other pairs are hard non-edges; two stubs cannot be together since then
their anchors would share a cluster although the pair of the anchors is an unedited hard non-edge). The comparator pad pairs inside `R_c` are plain edges there. All tunnels at c are unedited, so every anchor x is in a cluster `{x,stub}` (Lemma 2 of wave 6), every lock guest is with its z,
and by Prop. 4 (=>) with *all* tunnels uncut the position of c lies in no group, a contradiction (every position is forbidden by a lock or a permanent lock). QED. (VERIFIED: `w10_ex1.py`, l = 0 UNSAT, Sec. 6.5.)

**Theorem S (soundness; PROVED).** If `OPT(I_tile) <= |H| + k^2` then Grid Tiling (<=) is solvable. More precisely `ex(S) >= k^2` for every solution, and a solution with `ex = k^2` consists of one cut tunnel per cell
and exactly one edit per member and nothing else.
Proof. Let S be a solution, `ex = ex(S) <= k^2 = l`; p is large enough for Lemma 5, so no edit lies on a pad pair or on a pair inside a super-vertex; in particular the pad pair of every comparator pad member is unedited.
Counting. `|S| = sum_c |S cap Pi_c| - q + |S_rest|`, `S_rest` = edits outside all `Pi_c` (a pair in two `Pi` is a cross-anchor pair, counted once less). Every member of `H_cmp` has an edit (decomposition lemma),
and since its pad pair (inside some `R_c`) is not edited, its edit lies in `S_rest`: `|S_rest| >= |H_cmp|`. Let `h_c` = number of edits in `S cap Pi_c` that are on pairs of no member. Then `|S cap Pi_c| >= |H_c| + max(1, h_c)` (Lemma 6).
Summing: `ex >= sum_c max(1,h_c) - q`. Let `Z = {c : h_c = 0}`, `h_in` = hard edits inside a single `Pi_c`; then `sum_c h_c = h_in + 2q` and `|chains not in Z| <= sum h_c`.
So `ex >= |Z| + h_in + q` and `|Z| >= 2k^2 - h_in - 2q`; hence `ex >= 2k^2 - q`. Also `ex >= q` (cross-anchor pairs are hard). So `ex >= k^2`.
If `ex = k^2`: `q = k^2`, `|Z| = h_in = 0`, every chain has `h_c >= 1` and `sum h_c = 2q = 2k^2` so `h_c = 1`: the only non-member edits are `k^2` cross-anchor pairs, one at each chain. Such a pair is a tunnel edge: if an edited cross-anchor pair `(x_A(g), x_B(g'))` were a non-edge with `g' != sigma(g)`,
then `x_A(g)` is also tied to `x_B(sigma(g))` by its (unedited) tunnel, so two anchors of `B` would be together, but their pair is an unedited hard non-edge (`h_B = 1`): contradiction. So S = (k^2 deleted tunnels, one per chain, hence one per cell since tunnels join `A_ij` with `B_ij`) + (one edit per member).
Then S is a ZE partition of `I_tile - C`, and Prop. 4 (=>) gives a tiling. QED.

**Theorem C (completeness; PROVED).** If GT(<=) has a solution `(s_ij)`, then `OPT(I_tile) = |H| + k^2`. Proof: Prop. 4 (<=) gives a ZE partition of `I_tile - C`; its cost in `I_tile` is `|H| + |C| = |H| + k^2`; `ex >= k^2` by Theorem S. QED.

**Theorem H.** `(k,n,S) -> (I_tile, H, l = k^2)` is a polynomial reduction with `l + tw(I(H)) = O(k^6)` (and `l + c + pw = O(k^6)`), so CEaMP parameterised by `l + tw(I(H))` is W[1]-hard (GT(<=) is W[1]-hard in k). (Under ETH the same reduction gives no `f(l+w) n^{o((l+w)^{1/6})}` algorithm, since GT(<=) has no `f(k) n^{o(k)}` algorithm; this loss is due to p, w = poly(k); not needed for the dichotomy.)
Pair-maximality (no all-hard induced P3): VERIFIED-BY-COMPUTER for k=1, n=2, p=8 (`w10_pm.py`: 0 all-hard induced P3); by construction hard edges are glue cliques, thin hard edges (complete bipartite between super-vertices, so triangles) and tunnels (anchor has one hard edge).
Otherwise Lemma MAX of the paper turns it into a pair-maximal instance with `tw` increased by at most `l`.

## 6. Computer verification (all scripts in `code/`)

### 6.1 Comparator semantics (OR-switch guest), thin and thick
`w10_t1t.py` (thin) and `w10_t1.py` (thick, p = 8, n = 515, |H| = 3905): two chains m = 10, separator group u = 2 (positions 4,5), clause `(a <= 3) or (b >= 6)`. For all 100 pairs (a,b) the ZE feasibility equals
`a,b notin {4,5}` and `(a<=3 or b>=6)`; 0 mismatches thin, 0 mismatches thick (the thick run takes 826 s).
### 6.2 Whole reduction, tunnel-release semantics (`w10_red.py`, `w10_batch.py`)
For random Grid Tiling instances and **every** choice of one tile per cell (and for a sample of choices where the number of choices is large), ZE feasibility of `I_tile - C` equals "the choice is a tiling":
| k | n | instances | tests | mismatches |
|---|---|---|---|---|
| 2 | 2 | 5 seeds (`w10_red.py 2 2 s`) + batch seeds 100..159 | 36+8+36+72+27 exhaustive + 1,787 sampled | 0 |
| 2 | 3 | seeds 1,2 exhaustive (840, 270 choices) + batch 200..214 (833) | 1,943 | 0 |
| 3 | 2 | batch seeds 300..311 (18 chains, 608 sampled choices, 10 instances with a solution) | 608 | 0 |
Thick (p = 8, k = 2, n = 2, 2,916 vertices, 23,176 members, 45 s per solve): 4 choices (2 tilings, 2 non-tilings): all agree (`w10_thick_t.py`).
### 6.3 All cut sets of size <= k^2 (`w10_cuts.py`)
k = n = 2, 10 tunnels: all 386 cut sets of size <= 4: ZE(I - C) feasible for exactly the 15 sets that are tilings (one tunnel per cell, valid); 0 mismatches. (So no cut set of size < k^2 works and no non-tiling set of size k^2 works, in the TR model.)
### 6.4 Chain and profile certificates
`w10_chain.py` (explicit partitions ZE for m = 7..60, equal to the SAT enumeration for m = 7,8,9,10,13; transfer matrix); `w10_prof.py` (1,836 cases).
### 6.5 EXACT optimum on the real thick instances (all pairs editable) -- `w10_exact.py`
Lazy-triple SAT with a cardinality bound on the excess (hard-pair edits + member excess), p = 8. Results (UNSAT is a proof):
* k = 1 (n = 1 plus one extra end group, m = 8), chains A,B with one tile and its tunnel, |V| = 428, |H| = 3258: `OPT > |H|` (l = 0 UNSAT, 4 s), `OPT <= |H|+1` (l = 1 SAT, certified cost 3259, 116 s). `w10_ex1.py 8 yes`.
* Same without tile and tunnel (no sharing possible), |V| = 424, |H| = 3200: l = 0 UNSAT and **l = 1 UNSAT** (OPT >= |H|+2): each chain needs its own deviation (`w10_ex1.py 8 no`).
* Row gadget (two A-chains n=2 with comparator, anchors tunnelled to private stubs, l=2, |V|=725, |H|=5707; `w10_ex2.py`): NOT completed -- both runs (YES and NO instance) were still in clause-encoding setup (python seqcounter for 266k pair tokens) after 30 min of throttled CPU and were killed. No result is claimed for it; the comparator is covered by 6.1-6.3 and by the proof of Theorem S.

## 7. How the earlier obstructions were handled (Wave 8 Sec. 7 "host/parity/co-hosting bookkeeping")
* *Thin comparators do not exist* (pad regress, wave 8 Sec. 6): all comparator guests are thick; pads use perfect matchings of `K_p`, so parity is automatic (p even).
* *Co-hosting*: two guests in one host cluster are mutually hard non-adjacent. Hosts are `s_{4u+2}` only; two such hosts are never in one cluster (they are 4 spine steps apart), so permanent guests never meet. The separator host serves two comparator guests whose
  position sets (`a <= 2u-1` resp. `a >= 2u+2`) are disjoint. The slot hosts carry one lock guest each and no comparator. This replaces the failed "one guest per odd host cluster forbids both M_{h-1} and M_h" bookkeeping: a group is the pair `{2u,2u+1}` and the forbidden set is the whole group.
* *Parity of the threshold*: the comparator uses separator groups (two positions that are never used by a tile), so `a <= 2u-1` and `a >= 2u+2` partition the usable positions; the old need for "two centre members per side" disappears.
* *Equality vs inequality*: the OR-switch gives implications; Grid Tiling with <= needs only monotone constraints, so no equality comparator is needed.
* *Matching of tiles between the A and B orders* is done by hard tunnels (far hard pairs are free for `tw(I(H))`), exactly as in wave 6 Sec. 3.1.

## 8. Status table

| item | status |
|---|---|
| Lemma 1 chain: explicit ZE partitions P_a (all m) | PROVED (local finite check, by computer, m <= 60 + translation invariance) |
| Lemma 1(b) no other ZE partition: transfer matrix x^4(x-1)^2, N_m = m | PROVED by computer (finite certificate) |
| Lemma 2 profile classes / allowed sets | PROVED from Lemma 1(a); VERIFIED m = 8..60 |
| Lemma 3 blow-up (thick ZE = thin-guest ZE) | PROVED; VERIFIED thick vs thin (6.1, 6.2) |
| Prop. 4 tunnel-release semantics / comparators | PROVED; VERIFIED (6.1-6.3) |
| Lemma 5 no split (p >= l + 5) | PROVED |
| Lemma 6 region bound | PROVED |
| Theorem S soundness (ex >= k^2, structure at equality) | PROVED; exact tests (6.5) agree on small cases |
| Theorem C completeness | PROVED |
| Theorem H W[1]-hardness in l + tw(I(H)), l + c + pw | PROVED (modulo the finite certificates above) |
| exact OPT on full k=2 instance | NOT done (n ~ 3000 thick: lazy SAT infeasible here); replaced by 6.2/6.3/6.5 |
| hardness with c bounded (p constant) | OPEN (glue argument needs p > l) |
| hardness for c = 2 / soft-triangle-free | OPEN (excluded for soft-triangle-free by wave 9, Theorem A) |
