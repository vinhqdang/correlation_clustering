# Wave 8 (hardness prover C): W[1]-hardness of CEaMP in l + tw(I(H))? -- spines, chords, SAT/CEGIS gadget synthesis

STATUS OF THIS FILE: written incrementally; the section "Final status" at the end is authoritative.

Labels: PROVED = proof in this file; VERIFIED-BY-COMPUTER = exhaustive/exact finite check (script named), evidence only unless stated;
HEURISTIC; OPEN.  Code: `code/w8c_*.py` (library `w8c_lib.py`; pysat/CaDiCaL for exact ZE decisions, HiGHS for exact OPT).
Conventions as in wave6_hardness.md: soft pair = packed pair, hard pair = unpacked; members (a,b,c) with centre b; a solution is
zero-excess (ZE) if |S| = |H|; state L/R/M of a member = delete ab / delete bc / insert ac; "anchor" = vertex of a member with no hard edge,
"lock" = stub or tunnel on an anchor.

## 0. Summary (what is new)

1. (PROVED + VERIFIED) **One-hot selection chain OH_m exists at zero excess with bounded c and bounded width.**  An explicit family
   (n = 3m+2 vertices, |H| = 3m members, c = 4, every edge of G is packed, tw(I(H)) bounded) whose ZE solutions are *exactly* m solutions,
   one for each position a in {0..m-1}, with member-state word R^a M L^{m-1-a}.  This is the "persistent bit with a freely chosen flip
   position" that wave 6 could not find.  Sec. 3.  Certified for ALL m >= 7 by a finite-automaton argument (Sec. 3.3) and by exhaustive SAT
   enumeration for m <= 40.
2. (PROVED) **Lock lemma** (Sec. 1): a clean anchor can only *forbid* togetherness events "z ~ u" (u outside the soft neighbourhood of y);
   it cannot force or forbid separations.  This explains why every earlier attempt with the plain P3 spine (whose flip is a *separation*
   event: a singleton) failed, and why the *merge* event M of OH_m is the right polarity.
3. (VERIFIED-BY-COMPUTER, negative) Attachment obstruction: all SAT/CEGIS searches for (L) a lock gadget "anchor forbidding M at one position"
   and (C) a comparator between two chains, in the tested size classes, returned UNSAT or did not finish (Sec. 5).  Sec. 4 explains the
   mechanism ("pad regress": a vertex that joins chain clusters needs soft pairs to the whole cluster).
4. Conditional reduction (Sec. 6): OH_m + (L) + (C) would give W[1]-hardness; precise gadget specification.  NOT a proof of hardness.

(see "Final status" at the end for what is decided.)
