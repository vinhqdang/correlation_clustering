# Certificates

- `colab/`: the certificates and clusterings behind every checked bound and
  archived clustering reported in the manuscript (runs on the Colab
  machines), re-checked in `results/mpc/recheck.csv`.
- `../local/certs/`: the certificates of the runs made on local machines (the
  1800 s bound runs on every fourth PACE exact instance, and part of the star
  local search runs with the stopping rule of the KaPoCE branch-and-bound),
  re-checked in `results/mpc/recheck_local.csv`.

File names are `TAG_GRAPH_SEED.npz`, as for the result files in
`results/colab/` and `results/local/`; clusterings end in `.labels.npz`.  The
format is specified in `docs/certificate_format.md`.
