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
clustering, and writes `results/mpc/recheck.csv`.  It exits with a non-zero
status if any certificate is rejected.

## 4. Re-running the experiments

Each result file in `results/colab/` is produced by one command; the tag (the
first part of the file name) selects the script and the budget:

| tag | script | content | paper |
|---|---|---|---|
| `c2`, `c2x` | `colab_run_cert.py T SEED TAG GRAPH` | CertiFlip with checked certificate (SNAP 600 s, PACE exact 60 s) | Tables of Sections 7.2 and 7.3 |
| `lpack`, `ltri`, `lstar`, `lwarm` | `colab_run_lp.py T SEED TAG GRAPH` | bound baselines: sparse star packing, metric LP on P with triangle / star rows, bound without time limit | Sections 7.2, 7.3 |
| `s2h`, `s2h150`, `s2h60`, `s2p` | `colab_run_seq.py T SEED TAG GRAPH` | protocol v2 head-to-head PXMem vs KaPoCE (held-out SNAP, PACE heuristic) | Section 7.4 |
| `sccevo` (in `results/scc/`) | `run_scc.py T SEED evo GRAPH` | SCC on the complete encoding | Section 7.2 |

For example

```bash
RESULTS_DIR=out CC_ROOT=$PWD python experiments/colab_run_cert.py 600 0 c2 ca-GrQc
```

`experiments/colab_fleet.py` distributes these commands over Colab machines;
the job lists are in `results/colab/queue.json`.  The held-out graphs and the
analysis plan were fixed in `experiments/heldout.md` before the runs.

## 5. Tables and figures

```bash
python experiments/mpc_tables.py
cd paper_mpc && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## 6. Proofs

```bash
cd lean/CCProofs && lake build
```

builds the Lean 4 development of Appendix A.
