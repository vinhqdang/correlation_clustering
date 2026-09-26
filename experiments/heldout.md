# Held-out head-to-head: pre-registered plan

Written and committed before any run of protocol v2.  The results are
reported in the paper whatever they are.

## Instances

Twelve SNAP graphs that were never used while developing PXMem or CertiFlip
(no earlier run of either solver on them exists in `results/`).  Chosen to
span 10^3 to 10^6 vertices and several graph classes; the list is fixed:

| name | n | m | class |
|---|---:|---:|---|
| email-Eu-core | 1005 | 16064 | e-mail |
| facebook | 4039 | 88234 | social |
| wiki-Vote | 7115 | 100762 | voting |
| cit-HepTh | 27770 | 352285 | citation |
| cit-HepPh | 34546 | 420877 | citation |
| p2p-Gnutella31 | 62586 | 147892 | peer-to-peer |
| soc-Slashdot0902 | 82168 | 504230 | social |
| loc-Gowalla | 196591 | 950327 | location-based social |
| email-EuAll | 265214 | 364481 | e-mail |
| amazon0302 | 262111 | 899792 | co-purchase |
| web-NotreDame | 325729 | 1090108 | web |
| roadNet-PA | 1088092 | 1541898 | road |

(n, m: vertices and undirected edges after the instance definition of
Section 2, i.e. the positive graph without self-loops and duplicates.)

PACE 2021 heuristic track: the instances whose vertex count equals that of a
development graph are treated as development instances and excluded from the
held-out summary: heur167 (ca-GrQc), heur173 (ca-HepTh), heur175 (ca-HepPh),
heur176 (ca-AstroPh), heur178 (ca-CondMat), heur179 (email-Enron), heur180
(soc-Epinions), heur191 (com-DBLP), heur193 (com-Amazon), heur198
(com-Youtube).  The remaining 190 instances form the held-out PACE set.  (Seven
further PACE instances coincide with graphs of the list above; they are
held-out in both sets.)

## Protocol v2 (experiments/colab_run_seq.py, tags starting with `s2`)

* one machine per run; KaPoCE (PACE 2021 heuristic submission, seed patch and
  time limits from T as the only changes) and PXMem (default configuration,
  unchanged since the development runs) run one after the other, each pinned to
  the same logical CPU while the machine is otherwise idle;
* both read the same PACE-format text file inside their measured time;
* order of the two solvers random, fixed by (graph, seed);
* KaPoCE: internal limits 590/600·T and 90/600·T, SIGTERM at T;
* PXMem is scored at KaPoCE's elapsed wall-clock time from its full
  improvement history (it is also recorded at T);
* budgets T = 600 s (primary), 150 s and 60 s; seeds 0–9 on the twelve SNAP
  graphs; seed 0 on the PACE instances at 600 s;
* invalid KaPoCE output counts as a PXMem win (PACE rule) and is reported
  separately; a PXMem run without a solution at the cut counts as a loss.

## Analysis (fixed in advance)

* Unit of inference: the graph.  For each graph, the median over seeds of the
  relative difference (PXMem − KaPoCE)/KaPoCE.
* Primary test: two-sided Wilcoxon signed-rank test over the twelve graph
  medians at T = 600 s; also the sign test over graphs.
* Per graph: exact two-sided sign test over seeds with Holm correction over
  the twelve graphs (reported with its power limitation: with ten seeds only
  10/0 outcomes can reach significance), and the Hodges–Lehmann estimate of the
  paired difference with its exact confidence interval.
* PACE held-out set: one run per instance, so the instance is the unit;
  Wilcoxon and sign test over the 190 instances.
* Budget dependence: the same tests at 150 s and 60 s.
* Equivalence is not claimed from non-significance.
