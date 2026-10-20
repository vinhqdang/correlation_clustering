# Hostile referee report, round 1 (full read of paper/main.pdf, 51 pp.)

Scope: main.tex + all 15 section files + refs.bib, all proofs read; 18 proofs re-derived line by line (list in section 6);
brute-force scripts in `code/r1_*.py` (section 7). `paper/` untouched, no git. Locations are `file:line` of the .tex sources.

## 0. Verdict in five lines

* Mathematically the paper is in better shape than its form: every proof I re-derived is correct up to small gaps (listed in section 1); no explicit instance/claim failed a brute-force check.
* As a journal paper it is not acceptable in this form: the headline results (CSP view, `7^{w+1} n^{2l}` DP, ETH bound, Theorem K, `n^l` algorithm) occupy about 12 of 51 pages, and about 25 pages are partial, heuristic or negative material with five status tables, "the notes", "the audit", "not re-run for this paper".
* The central question (FPT in `l+w`) is open; the abstract nevertheless lists 12 items, many of them technical lemmas that nobody can use downstream.
* Novelty is thin where it is proved (the CSP/DP is immediate, the ETH bound is a corollary of the LPS23 reduction, Theorem K is the BFK18 rule applied to components) and the interesting part (T1/T2 locality, FAR obstruction, cluster-identity normal form, n^l) is not packaged as a theorem.
* Fix list (section 8) is prioritised; the restructured ~30-page outline is in section 5.

---

## 1. Correctness (statement vs. proof, hypotheses, constants)

### 1.1 Proof gaps / statements stronger than their proofs (all fixable, none fatal)

| # | Location | Problem | Fix |
|---|---|---|---|
| C1 | tunnels.tex:57-61, Lemma `tmaj` ("every core has at most two anchors") | Proof ends with "A binary subtuple of t involves at most two anchors". False as written: coordinates are *members*, and a pair of members can carry up to 4 anchors (a member has up to 2 end-anchors). What the Baker-Pixley argument actually gives is: if S is a minimal non-face then S = (anchors through two members), so `|S| <= 4`. To get 2 one needs the missing remark: a member with two end-anchors is never a face (both pairs at its centre would be apart = 2 edits in a ZE solution; the centre cannot be a singleton either), so such a pair is itself a core of size 2; a core of size >= 3 therefore has at most one anchor per member, and then a pair of members carries at most 2 anchors. Also say what happens when the forced partial vector t is inconsistent. | add the remark, or state bound 4 |
| C2 | reduction.tex:30-35, Prop `two` ("decides whether T is tight and whether it is reducible") | Proof covers only reducibility *with a tight S_T* (3-valued CSP). For a non-tight component the optimal S_T has > p_T edits and is not a 3^p assignment. Statement overclaims. | "for a tight component T, decides whether T is reducible" (this is all Thm `three` uses) |
| C3 | reduction.tex:46-52, Cor `three`(b) | Proof is one line ("the packing produced by the last reduction pass is already pair-maximal"). Missing: (i) Rule 6.2 preserves pair-maximality of H\F (true: after the rule no induced P3 of G+S_F has exactly two vertices in W_F because of the R-condition, and a P3 with <=1 vertex in W_F has the same pairs/ownership as in G); (ii) tw(I) does not increase when components are deleted; (iii) why the loop ends. | add a lemma with this argument (I verified it) |
| C4 | locality.tex:23-32, Thm `T1T2` ("exactly one of the following holds") | The witness z of Lemma `half`(b) is not unique, so "type" is not a property of the pair ab, only of (ab,z). Cor `softfree`(ii) uses it correctly (works with the chosen z) but Prop `cex` (a),(b) say "type (T2) and not (T1)" for a pair. | define: ab is of type T1 if *some* witness has both steps soft, else T2; (a),(b) then hold since Q(2)∩Q(5)=∅ |
| C5 | structure.tex:92-94, Remark tightness of shrink threshold | Tightness proved only for l=0 (C4, |X|=1). Abstract (v) says "a shrink rule with a tight threshold". I found a clean l=1 example: C4 with vertex 2 replaced by a clique X of two twins (n=5, H=(0,1,3)): OPT(G)=3, OPT(G-x)=2, so (G,H,1) is NO and (G-x,H,1) YES with |X|=l+1=2 (`code/r1_shrink.py`). For l=2 (|X|=3) the same family does *not* separate (OPT 3 vs 3), so general tightness is unproved. | state tightness for l in {0,1} only, or find a family |
| C6 | dproute.tex:88-89, Remark `rad2` | Lemma `half`(b) + NF(alpha) *prove* that inserted pairs are at distance exactly 2 in G[L_B]; the text then keeps rho=2D+l "because it was not part of the audited statements". This is not a mathematical reason. With rho=2 the bound `Z_*^{rho_* l}` with rho_*=2D_*+l=O(w+l) becomes `Z_*^{2l}`, i.e. Thm `nl` improves from (w+l)^{O(l(w+l))} to (w+l)^{O(l)} for free. Cor `cand`(ii),(iii) and Thm `nl` should be stated with rho=2. | do it |
| C7 | reduction.tex:12,17 vs dproute.tex:63 | Constant clash: Prop `one`(b) uses `Z* = max{Z,(l+1)(2D+2l+1)}`, but an exceptional cluster has size <= 2D+l+1 (NF(alpha)); the proof of (b) itself says "of size <= 2D+l+1" and then writes 2D+2l+1. | use (l+1)(2D+l+1) |
| C8 | excess.tex:1, intro.tex:13, main.tex abstract (ii), open.tex:4 | "XP in l, FPT in w" is not a parameterized-complexity statement: l is part of the input and the exponent depends on it, so the algorithm is *not* FPT in w; and Prop `w0` (limits.tex:15-20) says w alone is NP-hard. | say "XP in l with single-exponential dependence on w (`7^{w+1} n^{2l+O(1)}`)" |
| C9 | excess.tex:12-14 + abstract (ii) | "domain size 7 cannot be lowered in the natural encoding": only 3 is excluded (examples with 2 and 3 edits on one member, verified by me). Nothing on 4,5,6; "natural encoding" undefined. Also an encoding with domain 3 plus at most l 'multi-edit' members is possible. | delete the sentence from the abstract; keep the remark as "domain {L,R,M} does not suffice" |
| C10 | dproute.tex:148-150 + Thm `omega`/`P`, Prop `cand`, FPT claims in abstract (ix),(x) | Parameters (D, omega, tw(Gamma^l)) are those of a *pair-maximal* H', not of the input H, and a pair-maximal extension is not unique (greedy order), so `omega` is not a function of (G,H). Only Thm `nl` is stated for general H. Abstract (ix)/(x) say "for a pair-maximal H" but the title-level claim "FPT in (l,D,omega)" reads as a parameterization of the input. | state: "for every pair-maximal extension H' (e.g. the greedy one), time f(l',D',omega')n"; list in a table which theorem is for input H vs. H' |
| C11 | tunnels.tex:97, abstract (xii) | "a conjecture ... (majority closure for c<=2) is refuted": it refutes only the *uniform* operation m(all distinct)=M. Sigma={LLLL,LMLM,MLML,RRRR} (verified, `code/r1_maj.py`) may still be closed under another (non-symmetric) majority, e.g. choosing m(L,M,R)=L, m(M,L,R)=M gives LMLM. Lemma `tmaj` needs *some* majority operation. The conjecture is also one "of the working notes", i.e. unpublished. | say so; or drop the conjecture |
| C12 | locality.tex:15-18 vs abstract (viii) "(pieces) hence lie in one bag" | Prop `bag`(a) puts each *pair* (a,b) in a bag, not the pieces A, B; Cor `softfree` only gives A∪B inside a ball of radius 4 of W. | "each cross pair lies in a bag" |
| C13 | structure.tex:97, intro/abstract "size bound |R|<=(l+1)|V(H)|" | holds after exhausting the rules, for the pair-maximal H' (and `|V(H')|` is not bounded by l); abstract (v) omits both qualifiers. | add |
| C14 | intro.tex:17 | "examples with n=6,7,8 show that merges of pieces joined by a soft pair can be of type T2". Example (c), n=6, is *intra-piece* (locality.tex:89, 97). Abstract (viii) correctly says n=7 and n=8. | n=7,8 |
| C15 | pointers.tex:100-106 + abstract (x) "strictly stronger on a family" | Prop `Pstrict` needs in-degrees in {1,2} and l<=4; abstract has neither condition. | add or move |
| C16 | limits.tex:3-9, Thm `eth` | Proof relies on gadget sizes "read in the audit, not recomputed" and on "version-3 arXiv text". I re-read LPS23 (Sec. 4.1): cliques of size 5 (4m_i per variable), Q-cliques 1,4,14,1 and T-cliques 34,46,34 per clause => |V|=194m exactly as you state, so the claim is right. But it can be proved in two lines without gadget accounting: |V|=O(N+m) from LPS23's construction plus `3|H| <= 23|V|` (each vertex in <=23 members) gives |H|=O(m) and tw(I(H))<=|H|-1. Your "crude 3000m" is unnecessary. Also note the theorem is then equivalent to "no 2^{o(n)} algorithm for l=0" because tw<=|H|<=n: the tw statement carries no extra information (Remark (2) says this only implicitly). | rewrite proof; present as Corollary of LPS23 |
| C17 | locality.tex:15 wording / Lemma `near` | OK as stated (hypothesis "z,b in one cluster" is global to the lemma); but the pointers.tex:24 remark repeats the hypothesis "essential" -- consistent. No fix, listed because I checked it. | - |

### 1.2 Hypotheses used inconsistently across sections

* Sec. 5 (structure), 8 (locality), 9 (dproute), 10, 11: pair-maximal only; each section re-explains "otherwise apply Lemma MAX and read the primed instance" (structure.tex:97, packing.tex:45, locality.tex:5, dproute.tex:5, pointers.tex:5, reduction.tex:5) in slightly different words. Prelim Remark `max` (prelim.tex:17) also does it. One convention paragraph in Sec. 2 ("all results of Sec. X-Y are stated for a pair-maximal H; for input H use H' = H ∪ M, parameters w'<=w+|M|, c'<=c+|M|, D'=2c', l'=l-|M|") would do.
* Only Thm `nl` (dproute.tex:98-104) states primed constants (D_*=D+2l, w+l). Thms `omega`, `P`, Cor `cand`, Prop `candlist`, Prop `one`, Lemma `NF` do not; Remark `omegastate`(3) (dproute.tex:144) mentions it in passing.
* Sec. 6 (packing) is for arbitrary H except Lemma `HE` (pair-maximal); Thm `three` is "any packing" but Prop `one` needs pair-maximality; Sec. 11 says "none of the lemmas of this subsection and the next needs pair-maximality" (tunnels.tex:7) while Lemma `tfan` (tunnels.tex:65) assumes it. Make a table "lemma | needs pair-maximal | needs YES | needs optimal".
* Lemma SM, HUB, F1 etc. need "YES" for `u<=l`; some statements (locality.tex:53 `|A||B|<=u<=l`) silently assume it.
* `D`: (i) `D=D(H)=max|Q(v)|=2c` (packing.tex:4, Lemma `deg`), (ii) "any D >= max|Q(v)|" (SM, EX, NF, pointers), (iii) Cor `softfree`(iii) needs the equality D=2c to write `D(D-2)=4c(c-1)`, (iv) `D_*=D+2l`. Use `D` only for the number and `D(H)` for the function.

### 1.3 Constants that differ between sections

* exceptional cluster size: `1+2D+u <= 2D+l+1` (dproute.tex:63) vs `2D+2l+1` (reduction.tex:12,17) [C7].
* rho: `2D+l` (Cor `cand`, Thm `nl`) vs `2` (Remark `rad2`) [C6].
* `B=2D+3l`, `Z=B+D+1=3D+3l+1` consistent everywhere I checked (dproute.tex:56, tunnels.tex:66, reduction.tex:5, Prop `expander` margin "d_h<=4=2D+3l at l=0" ok). Thm `nl` uses `Z_*=3D_*+3l+1` with `D_*=D+2l`, i.e. `3D+9l+1`: consistent but the reader has to compute it; the abstract only says `f(l,w)`.
* Thm `K` (components.tex:5): `1.62^{(C(s,2)+1)l}`, kernel `(s^2-s+2)l`; checked: `4t_s+2 <= 2C(s,2)+2 = s^2-s+2` and s=5 gives 11l, 22l. OK. related.tex:2 quotes BFK18's `O(1.62^{(2t+1)l})`, same shape.
* ETH theorem constants: `|V|=194m` checked; "about 3000m" for |H| is just a crude bound, `1487m` follows from `3|H|<=23|V|`.
* Running times quoted in abstract/intro/open vs theorems: `7^{w+1}n^{2l+O(1)}` ok; `n^l f(l,w) 7^{w+l+1}` ok; `3^{w+1}n^{O(1)}` ok; `(Z(omega+1))^{O(omega)}n` ok.

### 1.4 Notation clashes (undefined or doubly defined)

* `c`: max members through a vertex (prelim.tex:8), number of edits of member `P_i` called `c_i` (dproute.tex:12) **while** `c_i` is also the third vertex of `P_i=(a_i,b_i,c_i)` (prelim.tex:4, expander, N*(F)), and `c(sigma)` the weight of a soft pair (pointers.tex:10).
* `W`: soft graph (packing.tex:109), `W_j` vertex set of a component (components.tex:2), `W_F` (reduction.tex:5), total weight `W=sum c(sigma)` (pointers.tex:35, Cor `LB`).
* `K`: kept hard edges (packing.tex:83), hard component (dproute), `K_v`, `K_v^+` cliques (dproute:152, reduction:59), `K_9`, `K_4`, Theorem K, `K_{m+1}`.
* `F`: family hard over H (packing.tex:4), deleted hard edges `F=U∩E` (packing.tex:83), cubic graph F (dproute.tex:149, reduction:59), tight-component items (F1)-(F6). **Two lemmas are titled "F1"**: dproute.tex:47 (clusters connected) and pointers.tex:42 (static pruning).
* `R`: `V\V(H)` (structure.tex:4), tunnel state `R` (tunnels.tex:7), R-condition, `R'`. `(R1),(R2)` = reducibility conditions (components.tex:38) **and** `(R1)-(R3)` = open items of locality.tex:114-117.
* `B`: budget `B=|H|+l-sum|x_i|` (excess.tex:16), cluster part `B_X` (structure), `B=2D+3l` (dproute.tex:56), shared members in a bag (tunnels.tex:44).
* `D`: `D(uv)` (structure.tex:19), `D_u,D_v` (structure.tex:11), Lemma D (reduction.tex:22), `D_*`.
* `Delta`: `Delta_1,Delta_2` (structure.tex:12), `Delta=B+D` (dproute.tex:77), `Delta(z,b)` (pointers.tex:7; the text warns).
* `I`: inserted pairs (packing.tex:83, locality.tex:5) vs `I(H)`.
* `T`: types T1/T2, components `T_j` (reduction.tex:5) vs `C_j` for the *same objects* in components.tex:1 (and `C` = cluster); tree `T`; `T=I(H)` in proof of (F2).
* `X`: class (structure), bag `X_t` (locality.tex:57), set in Cor `HUB`, members `X_i=(u_i,b_i,u_{i+1})` (reduction.tex:59).
* `M`: maximal hard family (packing), state `M` (tunnels), `M(z,b)` (locality.tex:76), `M_i` members (FAR), `M` in the proof of (F6).
* `l'`: `l-|M|` (packing.tex:22), `k-sum tau_j` (components.tex:19), `l-sum(tau_j-p_j)` (reduction.tex:38); primes `w',c',D'` only for MAX.
* `Q`: `Q(u)` soft partners vs quotient `Q'` (reduction.tex:13). `Gamma`: `Gamma_u`, `Gamma_h(z,b)`, `Gamma'`, `Gamma^l`. `N`: `N_h`, `N[Y]`, `N*(F)`, `N` = number of edges.
* "soft/hard" vs "packed/unpacked": defined in three places (packing.tex:4, locality.tex:5, dproute.tex:5, tunnels.tex:7); `\pairs` used, `\tw` used; `nu`, `alpha` (HUB) never defined; "Grid Tiling", "Baker-Pixley", "Robbins' theorem" uncited; Bocker's `O(1.62^k)` algorithm and the Chen-Meng `2k` kernel are used (components.tex:58,64) but **absent from refs.bib** and the text says "we have not opened the original sources" (components.tex:16).

---

## 2. Claims stronger than what is proved (abstract / intro / conclusion)

1. Abstract (ii) "domain size 7 cannot be lowered in the natural encoding" -- [C9].
2. Abstract (ii)/section title "FPT in w" -- [C8]; contradicts Prop `w0`.
3. Abstract (iii) "no `2^{o(tw(I(H)))}`": true but is just "no `2^{o(n)}` at l=0" in disguise and rests on unreproved gadget accounting -- [C16]. Say "ETH-hard in the exponent" and attribute the reduction.
4. Abstract (iv) "tw(G) and tw(I(H)) are incomparable": Prop `star`(b) (limits.tex:23,26) has tw(I(H))=0, tw(G)=k but with l = Theta(k^2) (H is far from maximal); it says nothing about the pair (l,w) being independent of tw(G). Only (a) (l+tw >= 2m-2) is meaningful.
5. Abstract (v) "a shrink rule with a tight threshold", "a size bound |R|<=(l+1)|V(H)|" -- [C5],[C13]. Also (v) is a side topic: nothing downstream uses Lemma R, Thm C1, Cor nf, Prop pointer (I grepped: Lemma `viab`, Prop `pointer`, Prop `B`, Lemma `exch` outside Sec. 5, Prop `c4`, Rem `c1` are never referenced elsewhere).
6. Abstract (viii) "hence lie in one bag" -- [C12].
7. Abstract (x) "strictly stronger on a family" -- [C15]; "an exact obstruction family" fine (Prop FAR is proved; I re-checked for m=3,4,5: 7/9/11 optima, all separate the chord).
8. Abstract (xi) "tight and reducible components decided in 3^{w+1}" -- [C2]; "so l=0 is fully reduced" is vacuous (Thm `l0` already solves l=0).
9. Abstract (xii) "refuted by an explicit instance" -- [C11]; "Computer evidence (not theorems)" inside an abstract is a symptom of the length problem.
10. Abstract closing paragraph and intro.tex:24 "What is not proved" repeat the same open items three times (also open.tex:4). The open question is stated honestly, which is good, but "A positive answer needs a bound g(l,w) on the number of candidate positions of the first endpoint" (open.tex:4) and "a hardness proof would have to use T2 or deletions" (locality.tex:121) are heuristics presented as necessities; the algorithm need not enumerate positions. Mark as "our route".
11. intro.tex:15 and abstract (vi): "proved ... in full" for Thm K but the load-bearing pieces (Lemma 6.3 of BFK18, Bocker, Chen-Meng) are cited; the abstract does flag this. Fine, but then Thm K is "BFK18 applied to components"; a referee will say so (see section 4).
12. Title: "When do local repairs ... combine into a cluster editing?" promises a characterization; the paper gives sufficient conditions plus an open question. Title page also reads "(working title; draft skeleton)" (main.tex:24) and the preamble has "% TODO: switch to the target template" (main.tex:7).

---

## 3. Structure and length

Page budget today (from the PDF): abstract 1; intro 2; prelim 1; CSP + XP DP (Sec. 3-4) 1.5; structure/twins/R/C1 5; Theorem K 2; packing counting lemmas 5; locality 5; dproute 8; pointers 5; reduction 5; tunnels 5; **limits (ETH, w0, star) 1**; open 1; related 0.5; bib 2.

* Headline does not stand as a theorem. There is no "Theorem 1" collecting: (a) ell=0 = CSP, scopes <= 3 variables = cliques of I(H), primal graph subgraph of I(H); (b) `7^{w+1}n^{2l+O(1)}`; (c) ETH; (d) failure of weaker conditions (Prop `w0`, vertex- vs pair-maximality, star trade-off); (e) Thm K; (f) `n^l f(l,w)7^{w+l+1}`. They sit in Sec. 3, 4, 6, 9, 13 with the ETH theorem on one page near the end.
* Roughly 25 pages are of the kind "what we tried": Sec. 7's last subsection (heuristic), Sec. 8 remarks `far`, `hardheur`, Sec. 9 `colour`, Sec. 10 (NFbad, Conjecture Lambda), Sec. 11 (TL, N*(F), 850-application evidence), Sec. 12 (all of it), plus five "Computer evidence" subsections and three status tables (dproute.tex:183-207, pointers.tex:145-163, reduction.tex:104-119; two of the tables cause overfull boxes of 16 and 22 pt).
* Process language must go: "the notes", "the audit", "audit's tests", "not re-run for this paper", "not audited again" (reduction.tex:93), "Semantic Scholar lists seven records (checked 2026-10-20)", "we did not re-derive", "read off the abstract of the arXiv version, 2026-10-20" (prelim.tex:8, related.tex:8-10, components.tex:16, tunnels.tex:94-96, packing.tex:119). A reference to `code/w5_packing_lemmas.py` (exists) should become a citation of an archived artifact.
* Move/cut list (exact):
  * Cut: Lemma `viab` (structure.tex:18-24, unused), Prop `B` + Remark (excess.tex:15-23, unused), Prop `pointer` (structure.tex:100-118, unused), Rem `c1` (structure.tex:126), Cor `LB` + BAD + its NO-test (pointers.tex:34-62; never used by a theorem), Prop `candlist` (pointers.tex:64-80; superseded by Cor `cand` and not used), Remark `colour` (dproute.tex:166), Remark `hardheur` (locality.tex:121), "Why simple hardness encodings fail" (packing.tex:121-130), Remark `far` (locality.tex:100), Conjecture `Lambda`, TL, `C2` (pointers.tex:140, reduction.tex:99, tunnels.tex:105), Problem `W2`, `twg`, `seth`, the whole "Related work" search log (related.tex:10).
  * Appendix (or companion/online section): Thm `omega` + DP + Prop `expander` (statement in main text as one proposition), Thm `P`/`Pstrict`, Prop `NFbad`, Prop `one`, N*(F) and its 3-page case analysis (reduction.tex:58-97; its (F6) proof was "repaired" and "not audited again"), the tunnel lemmas (keep only: lock lemma, memory lemma, TR-model proposition, chain-object problem; ~1.5 pages), twins/R/C1/N-forms (Sec. 5) as an appendix on kernel-type rules, SM/HUB/soft-graph counting lemmas (Sec. 6) as appendix.
  * Merge: Lemma `K3` (components.tex:47-54) and Thm `three` counting (reduction.tex:43) are the same argument written twice; state once as "Component counting lemma" (tight/reducible variants). Intro paragraph "What is proved" (11 lines, section-by-section) duplicates the abstract and should be a bullet list of 6 theorems.
* Abstract: currently ~1 page, 12 numbered items. Replace by <=200 words: CSP view + `3^{w+1}`, `7^{w+1}n^{2l}`, ETH, Thm K, `n^l` algorithm, "pair-maximal normal form and non-locality obstruction", question open.

---

## 4. Novelty / overclaim wording

* Honestly stated by the authors: twin lemma = critical-clique property (Guo09) -- but then "strict form" in the abstract/intro sounds like a contribution; Lemma `half` (half-degree, diameter 2) is a standard single-vertex-move fact in the Cluster Editing literature; cite instead of presenting as the "main tool". Theorem K is BFK18's rule applied to components (Lemmas K1-K3 are bookkeeping); the Corollary "I(H) perfect matching is FPT in l" is the visible consequence -- a referee will call it a corollary of BFK18 + the observation that components are induced subgraphs. State it that way.
* The CSP view and the DP are textbook (CSP of arity <=3 whose primal graph is a subgraph of I(H)). Presenting it as "What is proved" with a full DP proof is fine, but note the natural parameter is the treewidth of the primal graph (a subgraph of I(H), Obs. `two`), which is smaller and is what the algorithm really uses; the ETH bound is stated for I(H). Mention that tw(primal) <= tw(I(H)).
* Genuinely new and worth foregrounding: (i) the `n^l` improvement via the cluster-identity normal form (Lemma NF, Cor cand) -- with rho=2 [C6] this is clean; (ii) the T1/T2 dichotomy with the three small examples; (iii) Prop FAR (a deleted hard edge between members at distance m-1 in every optimum, w=1,c=2,l=1; verified); (iv) Prop `expander` (omega unbounded at l=w=0); (v) Prop `vertexmax` and C4 (pair- vs vertex-maximality).
* Wording to soften: "exact obstruction", "we prove ... five-vertex example shows that pair-maximality cannot be weakened" (true, verified) is fine; "FPT in w for each fixed l" (see C8); "clearly labelled heuristic" paragraphs inside Sec. 7/9/12 belong in the discussion, not in a results section.
* Dependence on a single cited lemma: Thm K and Thm `three` rest on Lemma 6.3 of BFK18. I read Rule 6.2 and Lemma 6.3 (arXiv 1512.04047): the rule is stated for an arbitrary *induced subgraph* H of G, matching your use; its proof is long and you do not reprove it. State clearly in the theorem statement "(uses [BFK18, Lemma 6.3])".

---

## 5. Suggested outline for a ~30-page journal paper

1. Introduction (3 p): problem, Main Theorem 1 (items a-f above, each with its section), related work (BFK18, LPS23, parameter choice), what is open; one status table (the only one).
2. Preliminaries (1.5 p): notation table (resolve all clashes of 1.4), soft/hard once, convention "pair-maximal H' / primed parameters" once, Lemma `decomp`, Lemma `deg`, Lemma `MAX` (needed everywhere).
3. Interaction CSP and the DP (3 p): Obs `one`,`two`, Thm `l0` (3^{w+1}), Thm `A` (7^{w+1}n^{2l}), Remark domain 7 (as "domain 3 does not suffice").
4. Limits (2.5 p): ETH corollary of LPS23 (proof via `3|H|<=23|V|`), Prop `w0`, Prop `star`(a), Prop `vertexmax` + C4 ("weaker maximality notions fail"), tw(I(H)) vs tw(G).
5. Bounded components (4 p): Thm K with one counting lemma (merge K3 and Thm `three`), Corollary perfect matching, scope remark. State explicitly what is cited.
6. Structure of optimal solutions for pair-maximal H (5 p): Lemma `half`, HE, H2, EX, `conn`, NF, Cor `cand` (with rho=2). Twin/R/C1/N-forms only as a half-page remark pointing to Appendix A.
7. The n^l algorithm (2 p): Thm `nl` with `Z_*^{2l}`; statement of Thm `omega` and Prop `expander` as a proposition ("the natural FPT route needs a bound on omega, which fails").
8. Non-local edits (3.5 p): Cor/Thm T1/T2 (type defined existentially), softfree corollary, Prop `bag`, Prop `cex` (a),(b), Prop FAR with proof, Prop `NFbad` (statement only).
9. Open problems (2 p): main question (one paragraph), chain object, c=2, base of exponent, treewidth of G.
Appendix A (not counted / online): twin lemma, R, C1, N1-N5, P1-P4; counting lemmas (SM, HUB, soft graph); Thm `omega` DP; Thm `P`; N*(F); tunnel section; evidence tables and code archive.

---

## 6. Proofs verified line by line (all correct unless noted)

1. Lemma `twin` (structure.tex:7-14): the two cost changes sum to exactly -2 (f(u,D)=f(v,D) for D avoiding u,v).
2. Lemma `exch` + Thm `C1` (structure.tex:42-56): arithmetic `x(x+2a-|C|)`, derived `a1>=x2+a2`, `a2>=x1+a1` contradiction -- correct.
3. Cor `nf` (N2)-(N5), shrink rule (N4) both directions, tightness remark at l=0 (C4 minus vertex 2 is the P3 itself) -- correct; tightness for l>=1 see [C5].
4. Prop `pointer` (P1)-(P4): P4 case split alpha>=q / alpha<q -- correct (but unused).
5. Lemma `HUB`(a)-(c) (packing.tex:64-70): in (c) |V(P)∩J'|<=2 because the middle vertex is adjacent to an element of J -- correct; Cor `HUB` correct.
6. Lemma `SM`(i)-(iii) (packing.tex:93-103): convexity step `r+u/r <= u+1` on [1,u], and the a<=D case `2|I| >= |C|(|C|-2D)` -- correct; `u<=l` needs YES.
7. Prop `soft` (treewidth <= 3w+2) -- correct (union of pairwise-intersecting subtrees).
8. Lemma `half`, Thm `T1T2`, Cor `softfree`(i)-(iii), Prop `bag` -- correct; [C4] on well-definedness of "type"; count `D(D-2)=4c(c-1)` correct for D=2c.
9. Lemma `EX`(a)-(c) incl. the del/ins bookkeeping (dproute.tex:38-43); Lemma `NF`(alpha),(beta),(gamma); Cor `cand` -- correct. Remark `rad2` is also correct and should be used [C6].
10. Thm `nl` (soundness via Remark `fixU`, completeness via Cor `cand`, time) -- correct.
11. Thm `omega` (both directions; DP introduce/forget/join formulae incl. `g1+g2-g(b)+x1x2`; connectivity => cluster closed when last bag vertex is forgotten) -- correct.
12. Prop `expander` (cubic F, three vertices per oriented edge, pair-maximality, d_h <= 4, minor) and Prop `Pstrict` (`|Delta(a_e,b_e)|=3+in(p)+in(q)>=5`) -- correct; I recomputed Delta on K4 with in-degrees (1,2,1,2): 6,6,5,6,7,6 as predicted.
13. Lemma `W1`, `W2` (<=2D per edited pair), Cor `LB`, `F1`, `DEL`, `BAD` (<=2c eps), Thm `P` (the dropped soft pairs have |Delta|>l hence are split) -- correct.
14. Prop `FAR` (i)-(iii) (pointers.tex:119-125): unicyclic count `f>=-1`, at most m clusters with f=-1, case analysis on {1,2m-1} -- correct; brute force m=3,4,5 confirms 7/9/11 optima.
15. Lemma `K1`-`K3`, Thm `K` (time, kernel arithmetic `4t_s+2<=s^2-s+2`) -- correct modulo BFK18 Lemma 6.3/6.4, Bocker, Chen-Meng (cited).
16. Lemma `D`, Thm `three` counting (|A|<=2l', |B|<=l-l', r<=l+l'), Prop `two` -- correct apart from [C2], [C3].
17. Lemma `tlock`, `tinactive`, `tsplice`, Prop `tTR` -- correct; Lemma `tmaj` [C1].
18. Prop `star`(a): OPT=2m-1 at j in {1,2} -- correct. Thm `eth` -- correct, see [C16].

---

## 7. Brute force (code/r1_*.py, all new, independent of the earlier scripts)

* `r1_examples.py` (all set partitions): domain-7 examples (2 and 3 edits on one member, unique optima) ok; Prop `vertexmax` (OPT=2, unique optimum {03,24}, H not pair-maximal) ok; Prop `c4` (three optima) ok; Prop `cex` (a) n=7 OPT 4, (b) n=8 OPT 6, (c) n=6 OPT 4, all unique optima, all valid and pair-maximal, as stated; Prop `NFbad` both n=10 instances (OPT 6, unique optimum, deleted edit 06 resp. 08, no locally bad member, I(H)-distance 2 resp. 3) ok; Prop FAR m=3,4,5 ok; `majority` instance Sigma={LLLL,LMLM,MLML,RRRR}, MMMM not in Sigma ok (`r1_maj.py`).
* `r1_gstar.py`: the 7-vertex flip gadget G* is valid, pair-maximal, OPT=5=|H|, exactly the five ZE solutions listed in tunnels.tex:94.
* `r1_random.py` (seeds 1,11-14; random graphs n=5-8, random maximal packings, all partitions): about 1.2 million checks, **0 failures**: Lemma R, Thm C1, (N2), (N3a/b), (P2)-(P4), Lemma half(a,b), T1T2, Lemma near, W1, W2, EX(a,b,c), DEL, SM(i),(iii), BAD, HUB(a),(c), decomposition identity, shrink (N4) for l=0..3.
* `r1_nf.py`: NF(alpha) (count, size, subset of L_B), connectivity, NF(beta) (no hard edge leaves, union of whole hard components, |C|<=2D) -- 0 failures. (BIG never occurs at n<=8, so NF(gamma) is only proved, not tested.)
* `r1_softfree.py`: Cor `softfree`(i),(ii) on 511 soft-free merges, 0 failures.
* `r1_shrink.py`: C4 with vertex 2 doubled: OPT(G)=3, OPT(G-x)=2 => shrink threshold l+1 fails at l=1 (new, supports C5).
* `r1_nstar.py`: N*(triangle) (n=18, |H|=8): valid, pair-maximal, no tight solution of the whole graph, tight on V\{y}, 95 pair-closed subfamilies, **none reducible** (consistent with (F6); the paper checked 90 with |F|<=6 plus the whole). `r1_expander.py`: K4 instance of Prop `expander`/`Pstrict` as predicted. Not run: N*(K4) (3^14 per subfamily).
* Not cheap / not done: Thm `omega` DP against brute force (the paper reports 2,778 agreements; I did not reimplement), Thm K against a cluster-editing solver.

---

## 8. Prioritised fix list (exact locations)

P0 (blocks submission)
1. main.tex:26-47 abstract: cut to <=200 words; remove items (vii)-(xii) as numbered claims.
2. Add a Main Theorem in the introduction (intro.tex:11-24 -> bullet list of Thms `l0`, `A`, `eth`, `K`, `nl`, `w0`/`vertexmax`), move ETH (limits.tex:3-9) to Sec. 4 with `3|H|<=23|V|` proof [C16].
3. Remove all process language (see section 3 list) and the three extra status tables; one table in the intro.
4. [C1] tunnels.tex:57-61; [C2] reduction.tex:30-35; [C3] reduction.tex:46-52 -- add the missing arguments.
5. [C6] dproute.tex:77-90, 98-104: state Cor `cand` and Thm `nl` with rho=2, delete "was not part of the audited statements".
6. [C8],[C9],[C10]: excess.tex:1, intro.tex:13, abstract (ii),(ix),(x), open.tex:4: wording of XP/FPT, domain 7, primed vs input parameters.
7. components.tex:16,58,64: add Bocker (1.62^k) and Chen-Meng (2k kernel) to refs.bib; drop "we have not opened".
8. Cut list of section 3 (unused lemmas, heuristics, conjectures) and move the rest to an appendix.

P1
9. [C7] reduction.tex:12,17 constant; [C5] structure.tex:92-94 and abstract (v) (use the l=1 example above); [C13],[C15],[C14] (intro.tex:17: n=7,8).
10. [C4] locality.tex:23-32 and 87-88: define type existentially. [C12] abstract (viii).
11. [C11] tunnels.tex:97-108, abstract (xii): "uniform majority".
12. Resolve notation clashes of 1.4; rename the second "F1" (pointers.tex:42 -> "static pruning", dproute.tex:47 -> `conn`); rename (R1),(R2) of locality.tex:114-117 to (O1)-(O3); `C_j` vs `T_j` for components (components.tex:1 vs reduction.tex:5); one definition of soft/hard.
13. Convention paragraph on pair-maximal/primed parameters in Sec. 2; one table lemma | pair-maximal | YES | optimal.
14. Abstract/intro: attribute Lemma `half` and the twin lemma, and Thm K = BFK18 on components (section 4).

P2 (typesetting, bibliography)
15. Overfull boxes (log): components.tex:12-13 (2.8pt), packing.tex:76-77 (1.5pt), locality.tex:88-89 (2.5pt, Prop `cex`(b) edge list: use `\allowbreak` or a table), dproute.tex:184-207 (16.5pt table), pointers.tex:9 (5.9pt display: split the definition of `Delta,Gamma_h`), pointers.tex:146-163 (22pt table), tunnels.tex:39-40 (8pt), related.tex:8-9 (3.3pt); underfull (badness 10000) tunnels.tex:94-95 (the `\sloppy` item with long edge lists).
16. main.tex:7 "TODO", main.tex:24 "(working title; draft skeleton)", title promises a characterization.
17. Theorem-like names `prop:one`, `thm:three`, `prop:five` are leftovers (labels only, harmless); counter shared by Theorem/Lemma/Proposition leads to numbers up to ~90 (log shows "Lemma 87"); consider per-section numbering.
18. refs.bib: LPS23 is cited as ACM TALG 2023 in the bib but the text says "version-3 arXiv text" (limits.tex:11 (4)); Korhonen entry carries a verification note; add arXiv/DOI only, no "metadata verified" notes in a bib file shipped with the paper.
19. Related work (related.tex:10): replace the search log by two sentences ("to our knowledge no prior work ...").
