# Analysis 1: the l = 0 case is a bounded-domain CSP on the intersection graph

Notation. H = {P_1..P_p}, P_i = (a_i, b_i, c_i) an induced path a-b-c with
edges ab, bc and non-edge ac; pairs(P_i) is its 3 vertex pairs. Modification-
disjoint means pairs(P_i) are pairwise disjoint. Intersection graph I(H): nodes
H, P_i ~ P_j iff V(P_i) and V(P_j) share a vertex.

## Observation 1 (proved)  [l = 0 is a CSP]
Every solution with l = 0 uses exactly one edit inside each pairs(P_i) and no
other edit (each P_i needs >= 1 edit among its own pairs, these pair sets are
disjoint, and the budget is |H|). So a solution is a vector x in {1,2,3}^p
(which pair of P_i is edited). G + S(x) is a cluster graph iff no triple of
vertices {u,v,w} induces a P3. For a triple, the status of each of its 3 pairs
is either fixed (the pair is in no packed P3) or a function of one x_i. Hence
the problem is a CSP with p variables of domain size 3 and one constraint per
triple, whose scope is the set of owners of the triple's packed pairs.

## Observation 2 (proved)  [scopes are cliques of I(H)]
Two pairs of one triple share a vertex. If uv in P and vw in P' then v is in
V(P) and V(P'), so P ~ P' in I(H). Hence every constraint scope is a clique of
I(H) (at most 3 variables). Consequently a tree decomposition of I(H) of width
w yields a dynamic program over bags in time 3^(w+1) * n^O(1) for l = 0.

Checked by computer: `code/check_csp_view.py` (486 random instances on 4..7
vertices) confirms Observation 2 on every triple and the equivalence "OPT = |H|
iff some choice vector gives a cluster graph".

## What this changes about the proposal
The proposal says that two disjoint P3s can still interact through the
transitivity of clustering, so a graph on P3s sharing vertices is not enough.
For l = 0 that is false by Observation 2: it is enough, and w = tw(I(H)) gives
a 3^O(w) algorithm. The interaction the proposal worries about is real but it
comes from the excess edits (l > 0):
* an excess edit may be on a pair outside all packed P3s, between vertices that
  lie in no common P3 (it lies in the scope of n - 2 triples);
* a packed P3 may receive 2 or 3 edits.
Guessing the unpacked excess edits costs n^(2l). That is the XP algorithm for
c = 2; it is also exactly where c = 1 needs the FPT argument of van Bevern et al.

Consistency checks against known results:
* c = 1: I(H) has no edges (w = 0), still non-trivial for l > 0. So the
  difficulty of l > 0 cannot be a function of w alone.
* c = 23, l = 0: I(H) has maximum degree <= 66 but unbounded treewidth. So
  Theorem 1 does not contradict Observation 2.

## Open problems that this analysis points to (none solved)
1. FPT in (l, tw(I(H))): needs a way to avoid guessing the l unpacked edits.
   Candidate: show that an excess edit on an unpacked pair can be restricted to
   pairs within bounded distance of V(H) in a graph built from I(H) plus G,
   by an exchange argument; or show W[1]-hardness in (l, tw), which would be
   the 'counterexample' part of the paper.
2. Is c = 2 FPT in l? (Open in Li et al. Observation 2 reduces it to handling
   the excess edits on a subcubic intersection graph.)
3. Lower bounds: is 3^o(w) impossible for l = 0 under ETH? Needs a reduction
   from a 3-valued CSP / 3-colouring on bounded-width instances to CEaMP.
4. Does a parameter coarser than tw(I(H)) suffice (e.g. one defined on the
   intersection graph restricted to vertices of G of degree > 2)?
   Real graphs: packings from star packings put thousands of P3s on a hub, so
   I(H) has huge cliques and tw(I(H)) is large; a useful parameter must treat
   hubs differently. This is a question for the experiments, after the theory.

## Next steps
* Read sections 3-5 of arXiv 1910.08517 (the c = 2 algorithm) to see exactly
  which step is XP and why.
* Citation search for follow-up work (Semantic Scholar, Google Scholar).
* Attempt open problem 1: first for c = 2, where I(H) is subcubic.
