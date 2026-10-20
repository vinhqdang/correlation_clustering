# Wave 8 (prover C): CEaMP in l + tw(I(H)) -- chain objects, lock lemma, status of the hardness attempt

Labels: PROVED = argument in this file; VERIFIED-BY-COMPUTER = exact finite check (script named; evidence for the instance sizes tested only);
HEURISTIC; OPEN. Code: `code/w8c_*.py` (library `w8c_lib.py`; pysat/CaDiCaL for exact ZE decisions, HiGHS for exact OPT).
Conventions as in wave6_hardness.md: soft pair = packed pair, hard pair = unpacked; member (a,b,c) has centre b; ZE (zero excess) = exactly one
edit per member and no other edit; member states L/R/M = delete ab / delete bc / insert ac; anchor = member vertex with no hard edge;
lock = forcing an anchor to be a singleton.

## 0. Verdict

The open question (W[1]-hardness vs FPT of CEaMP in l + tw(I(H)), or l + c) is NOT decided by this wave. What was obtained:

* the missing "chain object" EXISTS (one-hot selection chain, Sec. 3), with bounded c and bounded width;
* a polarity theorem (Lock Lemma, Sec. 1) that explains every earlier failure and says which attachments are possible;
* bounded 3-cores exist (Sec. 4), contradicting the wave-6 evidence that cores are large;
* a thick (blown-up) lock/permanent-guest construction that gives a genuine N-core chain verified for N=5 (Sec. 5);
* the comparator between two chains (the part needed for a Grid-Tiling / Multicolored-Clique reduction) is NOT built: OPEN.
  No impossibility theorem for it was found either, so nothing here points to FPT; it also does not give hardness.

## 1. Lock Lemma  (PROVED; numerically cross-checked, `w8c_lockslem.py`, 0 failures)

Let x be a clean anchor, an end of Q=(x,y,z) with centre y. Locking x (x a singleton) forces Q into state L (delete xy), hence y~z remains.
Consequences:
1. A lock forbids only togetherness events "z ~ u" for u outside the soft neighbourhood of y; it cannot force or forbid separations.
2. y must be closable (y's own cluster must be able to absorb z).
3. Moving y into C(z) is valid iff every pair (y,c), c in C(z), is soft.
So the plain P3 spine (whose flip event is a separation, a singleton) cannot be read by a lock; the flip must be a MERGE event.
This is why the Prop FAR / spine-with-chord idea of angle (A) fails in its naive form: the cut chord is a separation event.

## 2. Pad regress  (PROVED, short)

A guest vertex g joining a cluster must have soft pairs to every vertex of that cluster, so each joinable cluster needs soft pairs to all its
vertices; the resulting pad members use always-together intra-cluster pairs and are cost-neutral. A guest with k single covers needs
|W| = k (mod 2). This is the mechanism behind the earlier UNSAT results for thin gadgets (Sec. 6).

## 3. One-hot selection chain OH_m  (explicit; m ZE solutions)

Family `oh_instance(m)` in `code/w8c_oh.py`: n=3m+2 vertices, |H|=3m members, c=4, all edges packed (pure soft), bounded width.
* VERIFIED-BY-COMPUTER: ZE solutions are exactly m, member-state word R^a M L^{m-1-a}, a=0..m-1, for every m=7..16, 20, 30, 40
  (`w8c_oh.py`, `w8c_transfer.py`).
* VERIFIED-BY-COMPUTER (certificate): a finite-automaton / transfer-matrix certificate (`w8c_transfer2.py`) covers all m>=7 via the linear-recurrence
  argument. Treated as a proof sketch, not a written-out proof: label HEURISTIC-STRONG for general m.
This is the "persistent bit with freely chosen flip position" requested in the task; the selected position is where the merge M occurs.

## 4. Small cores  (VERIFIED-BY-COMPUTER)

CEGIS search (`w8c_core.py`, `w8c_per.py`): a 3-core exists with n=12, |H|=8 and c=4 (also c=3 variants): every 2 of 3 anchors can be locked
simultaneously by a ZE solution, all 3 cannot. Contradicts the wave-6 expectation that cores must be large. Exhaustive K=4 searches with n=14..16
(logs in /tmp/w8c) were unfinished.

## 5. Thick lock gadgets and the pair-chain PFG  (VERIFIED-BY-COMPUTER at one size; general statement HEURISTIC)

Blow-up lemma (argued, not fully written): replacing each vertex by a glued K_p proto-cluster and each member by p^2 Latin members gives ZE
solutions in bijection with the thin instance, and lets soft pair requirements of Sec. 2 be met by pad members. With permanent guests and lock
guests (`w8c_chainbuild.py`): chain m=9, p=8, n=247, |H|=1961, c=34: unlocked instance tight; all five lock groups simultaneously infeasible;
every "all but one group" feasible (M at positions 0,2,4,6,8), `w8c_chainbuild.py pairpfg`, SAT via `ze_solve_fast`.
So a persistent-flip chain with N=5 independently lockable positions exists, c=34 (bounded in terms of p only, not in m or N by design).
Generality in m follows from OH (Sec. 3), the blow-up lemma and the guest lemma, but is not written as a theorem.

## 6. Negative / unfinished searches

* Thin comparator between two OH chains (`w8c_cmp.py`): UNSAT or timeout in all tested size classes (VERIFIED-BY-COMPUTER, evidence only);
  consistent with pad regress.
* Thin lock gadgets on the plain spine (`w8c_det.py`): UNSAT, consistent with Lock Lemma.
* Periodic / nospine searches (`w8c_per.py`): found OH-type chains; no comparator.

## 7. Conditional reduction  (HEURISTIC / OPEN)

If a thick comparator existed (accept exactly the pairs (a_U, a_A) with [a_U>t]=[a_A>t] at even thresholds, built from comparator guests accepting
only X/Y clusters at even-index hosts, two centre-members per side for parity), then Grid Tiling (l=k^2, one value track per row/column and a tile
track per cell, tunnels matching pi, stub locks) would give W[1]-hardness in l + pathwidth, by the accounting "excess >= #tracks/2 forces the
tunnel-release form". Missing: (i) comparator verified at scale (two chains, n~250-400, design stalled on host/parity/co-hosting bookkeeping:
co-hosted guests are mutually hard non-adjacent, so one guest per odd host cluster must forbid both M_{h-1} and M_h), (ii) soundness proof.
None of this is proved.

## 8. Final status

* PROVED: Lock Lemma, pad-regress parity.
* VERIFIED-BY-COMPUTER: OH_m (m<=40, automaton certificate all m>=7), 3-cores, thick PFG N=5.
* HEURISTIC: blow-up and guest lemmas in general; conditional Grid-Tiling reduction.
* OPEN: comparator gadget, full reduction, soundness, and the FPT/W[1] dichotomy itself. No obstruction theorem against the comparator was found;
  next step is the multi-chain thick builder with comparator guests (`w8c_chainbuild.py` extended, verified by `ze_solve_fast`).
