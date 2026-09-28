# Reproducing the results of `paper_mpc/`

All numbers and tables of the manuscript are generated from the raw result
files in `results/` by `experiments/mpc_tables.py`; every lower bound comes
from an archived certificate in `results/certificates/colab/` and can be
re-checked without running any solver.

## 1. Environment

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements-lock.txt && pip install -e .
```

The runs used Google Colab machines (2 logical cores, Intel Xeon 2.20 GHz,
about 13 GB of memory), one run per machine at a time.  External baselines:

* KaPoCE, branch `heuristic_submission`, commit `64e2101` of
  https://github.com/kittobi1992/cluster_editing, built as in
  `experiments/kapoce/README.md`; `experiments/colab_run_seq.py` applies the
  only changes (seed and time limits from environment variables);
* SCC, commit `05d0849`, see `experiments/scc/README.md`.

## 2. Data

```bash
python experiments/fetch_data.py
```

downloads the 27 SNAP graphs and the 400 PACE 2021 instances (repository
commit `67f5df3`) into `data/raw/` and checks all 427 files against
`data/MANIFEST.sha256`.

## 3. Re-checking the certificates (no solver needed)

```bash
bash experiments/recheck_all.sh
```

re-parses every instance with the independent reader of
`experiments/check_certificate.py`, checks the instance hashes, every row and
the exact value of every certificate, recomputes the cost of every archived
clustering, and writes `results/mpc/recheck.csv` (Colab runs) and
`results/mpc/recheck_local.csv` (`results/local/certs/`).  The checker can
also be run without its compiled kernels (`NUMBA_DISABLE_JIT=1`), much more
slowly.  The certificate format is specified in `docs/certificate_format.md`.  It exits with a non-zero
status if any certificate is rejected.

## 4. Re-running the experiments

Each result file in `results/colab/` is produced by one command; the tag (the
first part of the file name) selects the script and the budget:

| tag | script | content | paper |
|---|---|---|---|
| `c2`, `c2x` | `colab_run_cert.py T SEED TAG GRAPH` | CertiFlip with a checked certificate (SNAP 600 s, PACE exact 60 s) | Sections 7.2, 7.3 |
| `lsstar` | `colab_run_lp.py T SEED lsstar GRAPH` | star local search (`experiments/sstar/sstar.cc`, compiled on first use with `g++ -O2`), runs until T; SNAP 600 s, PACE exact 60 s, seeds 0-4 | Tables 2, 6; Sections 7.2, 7.3 |
| `lsstar+k100` | same, tag `lsstar+k100` | as `lsstar` with at most 100 leaves per star (`SSTAR_KMAX=100`); com-Youtube, email-EuAll | Table 6 |
| `lsstar+kr` | same, tag `lsstar+kr` (or `run_local_sstar_kr.py`) | star local search with the stopping rule of the KaPoCE root bound, T = 60 s only as a cap; PACE exact, seeds 0-4 | Table 2 |
| `lsstar+t1200`, `lsstar+t3600`, `lsstar+k100+t3600` | same, T = 1200 / 3600 | equal-budget control for the block LPs, long runs on unconverged graphs | Section 7.3 |
| `lsslp` | same, tag `lsslp`, T = 1200 | star local search for T/2, then block LPs for T/2 | Table 6 |
| `lpack`, `lplp`, `lpack+eq` | same | ruin-and-recreate packing; packing then block LPs; packing alone with the budget of both | Table 6, Appendix B |
| `ltri`, `lstar`, `lwarm` | same | metric LP on P with triangle / star rows; bound with 1800 s and one block | Tables 2, 7 |
| `lkroot`, `lkstar` | same (modes `kroot`, `kstar`; needs the KaPoCE build of `experiments/kapoce/`) and `run_local_kroot.sh` | root bounds of the KaPoCE branch-and-bound on the input graph; its star packing exported and checked | Section 7.2, Table 6 |
| RAMA | `experiments/rama/` (GPU) | multicut dual on the support | Table 7 |
| `s2h`, `s2h150`, `s2h60`, `s2p` | `colab_run_seq.py T SEED TAG GRAPH` | protocol v2 head-to-head PXMem vs KaPoCE (held-out SNAP, PACE heuristic) | Appendix C |
| `sccevo` (in `results/scc/`) | `run_scc.py T SEED evo GRAPH` | SCC on the complete encoding | Section 7.2 |

The star local search stops on wall-clock time (except with `+kr`), so a
rerun with the same seed gives a similar but not identical packing; its
certificate is archived and re-checkable.  Machines: Colab VMs with two
logical cores (Intel Xeon 2.20 GHz on most runs, AMD EPYC on some; the CPU
model is in every result file, and the compiler in the files of the
revision runs), and a 4-core workstation for the rows marked in the paper.
Runs with the tag `+kr` that are not in `results/colab/` were made on a
4-core cloud container (`results/local/`).

For example

```bash
RESULTS_DIR=out CC_ROOT=$PWD python experiments/colab_run_cert.py 600 0 c2 ca-GrQc
```

`experiments/colab_fleet.py` distributes these commands over Colab machines;
the job lists are in `results/colab/queue.json`.  The held-out graphs and the
analysis plan were fixed in `experiments/heldout.md` before the runs.

## 5. Tables and figures

Two commands regenerate every table, number and figure of the manuscript.
They start from the archived result files and certificates in `results/`;
neither reruns a solver, except `make_figures.py`, which runs the solver code
on a small neighbourhood of ca-GrQc:

```bash
python experiments/mpc_tables.py all     # tables/*.tex, numbers.tex, two figures
python experiments/make_figures.py       # figures 2, 3 and 5
cd paper_mpc && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

| output | written by | from | runs behind it (Section 4) |
|---|---|---|---|
| Table 1 `tables/instances.tex` | `mpc_tables.py` (`instances`) | `results/mpc/instances.json`, `data/MANIFEST.sha256` | `experiments/instance_meta.py` |
| Table 2 `tables/pace_bounds.tex` | `mpc_tables.py` (`pace_exact`) | `results/colab`, `results/local`, `results/pace_exact.csv`, `data/pace2021_exact_kapoce_bounds.csv` | `lpack`, `c2x`, `lwarm`, `lsstar`, `lsstar+kr`, `ltri`, `lkroot`; triangle packings from `run_bench.py` |
| Table 3 `tables/pace_density.tex` | `mpc_tables.py` (`pace_exact`) | as Table 2 | as Table 2 |
| Table 4 `tables/pace_open.tex` | `mpc_tables.py` (`pace_exact`) | as Table 2 | as Table 2 |
| Table 5 `tables/pace_primal.tex` | `mpc_tables.py` (`pace_primal`) | `results/colab`, `results/scc`, `results/pace_exact.csv` | `c2x`, `sccevo`, `run_bench.py` |
| Table 6 `tables/snap_bounds.tex` | `mpc_tables.py` (`snap`) | `results/colab`, `results/snap*.csv`, `results/kapoce_root_snap_*.csv` | `c2`, `lpack`, `lplp`, `lsstar`, `lsstar+k100`, `lsslp`, `lkroot`, `lkstar`, all `lsstar` seeds and long runs |
| Table 7 `tables/snap_lp.tex` | `mpc_tables.py` (`snap`) | `results/colab`, `results/rama_snap_user.csv` | `ltri`, `lstar`, RAMA |
| Table 8 `tables/snap_sls.tex` | `mpc_tables.py` (`snap_sls`) | `results/colab` | `lsstar` seeds 0-4, `lsstar+t1200`, `lsslp`, `lsstar+t3600` |
| Table 9 (Lean coverage) | written in `appendix.tex` | `lean/CCProofs` | `lake build` |
| Table 10 `tables/snap_eqtime.tex` | `mpc_tables.py` (`snap`) | `results/colab` | `lpack+eq`, `lplp` |
| Tables 11-13 `tables/h2h_*.tex` | `mpc_tables.py` (`h2h`) | `results/colab` | `s2h`, `s2h150`, `s2h60` |
| Table 14 `tables/pace_heur2.tex` | `mpc_tables.py` (`pace_heur`) | `results/colab` | `s2p` |
| Figure 1 (overview) | `figures/overview.tex` (TikZ, static) | -- | -- |
| Figures 2, 3, 5 `fig_support`, `fig_gapmap`, `fig_crossover` | `make_figures.py` | solver code on a neighbourhood of ca-GrQc | computed when the script runs |
| Figure 4 `fig_sls_anytime.pdf` | `mpc_tables.py` (`anytime_figure`) | trajectories in `results/colab/lsstar+t3600_*.json` | `lsstar+t3600` |
| Figure 6 `profile_pace_heur.pdf` | `mpc_tables.py` (`pace_heur`) | `results/colab` | `s2p` |
| every number in the text | `mpc_tables.py` → `paper_mpc/numbers.tex` | as above, plus `results/mpc/recheck*.csv` | `recheck_all.sh` |

`mpc_tables.py` also asserts that every value accepted by the complete
re-check (Section 3) equals the value checked during its run.

## 6. Proofs

```bash
cd lean/CCProofs && lake build
```

builds the Lean 4 development of Appendix A.
