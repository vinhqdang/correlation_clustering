# Lower-bound certificate format (version 1)

This document specifies the certificates that `experiments/check_certificate.py`
accepts. It is written independently of the checker's code, so that other
solvers can write certificates without reading that code. Section 4 of the
paper gives the mathematics.

## 1. Instance semantics

A certificate refers to one raw instance file, exactly as it is distributed.

**PACE `.gr` files** (format `pace`):
- Lines that are empty or start with `c` are ignored.
- The header `p cep n m` gives `n`.
- Every other line `a b` is an edge between the 1-indexed vertices `a` and `b`. Vertex `a` has index `a-1`.

**SNAP edge lists** (format `snap`, optionally gzip-compressed, `.gz`):
- Lines that are empty or start with `#` or `%` are ignored. Commas count as whitespace.
- Columns 0 and 1 are the endpoint ids of an edge.
- In a signed file with sign column `s >= 0`, only the lines whose column `s` is a positive number are kept.
- The vertices are the distinct ids that occur in a kept line, including a line with equal endpoints. They are numbered `0..n-1` in increasing id order.

**Both formats**:
- Edges are undirected.
- Self-loops are dropped, and duplicate edges are merged.
- The result is the positive graph `G+ = (V, E+)` of a complete instance: every pair that is not an edge is a negative pair.

**Distance-two support**: `P` is the set of pairs `{u,v}` that are an edge or have a common neighbour. All other pairs are *far*.

## 2. File

The certificate is a NumPy `.npz` archive written with `allow_pickle=False`. It contains the following arrays; scalars are 0-dimensional arrays.

### Identity (mandatory)

| key | type | meaning |
|---|---|---|
| `n` | int | number of vertices |
| `instance_m` | int | number of edges after canonicalisation |
| `raw_sha256` | str | SHA-256 (hex) of the raw file's bytes, as stored (compressed if `.gz`) |
| `edge_sha256` | str | SHA-256 (hex) of the canonical edge list, see below |
| `instance_file` | str | path of the raw file relative to the data directory, e.g. `pace/exact/exact004.gr` (required in directory mode) |
| `instance_format` | str | `pace` or `snap` (optional; in directory mode it must agree with the instance table) |
| `sign_col` | int | sign column, `-1` if unsigned (optional; as above) |

**Edge hash.** Order the edges `(u, v)`, with `u < v` as 0-based indices, lexicographically by `(u, v)`. The hash is SHA-256 over:
- `n` as one little-endian signed 64-bit integer, followed by
- the pairs as consecutive little-endian signed 64-bit integers `u0 v0 u1 v1 ...`.

### Rows and multipliers

The rows are stored in compressed-sparse-row form, in the x-form `sum_p a_ip x_p >= b_i`. Here `x_p = 1` if the pair `p` is separated.

| key | type | meaning |
|---|---|---|
| `ptr` | int, length `R+1` | row `i` has entries `ptr[i]..ptr[i+1]-1`; `ptr[0]=0`, non-decreasing, `ptr[R]` = number of entries |
| `u`, `v` | int | the two vertices of each entry's pair (0-based; order irrelevant) |
| `val` | float, integral | coefficient `a_ip` of each entry |
| `b` | float, integral | right-hand side `b_i` of each row |
| `y` | float | multiplier `y_i >= 0` of each row, finite |
| `bound` | float | the value the solver claims (optional; reported only) |

### Constraints on the data

- Every array is one-dimensional, and none may be missing.
- `ptr`, `u` and `v` have an integer type.
- `val`, `b` and `y` are real: integer or floating point, never complex.
- `u`, `v` and `val` have the same length.
- `b` and `y` have one entry per row, and `ptr` has one more.
- Every value is finite.
- `0 <= u, v < n`.

- `|val| < 2^52` and `|b| < 2^52`, and both are integral.
- `y * 2^30 < 2^62`.
- No row contains the same pair twice.
- The absolute coefficients of each row sum to less than `2^61`.

## 3. Meaning and check

Let `c_p = +1` for an edge and `-1` for a negative pair of `P`. The checker proceeds in six steps.

1. It checks the identity fields against the file, as parsed by the rules above.
2. It checks that every pair of every row is in `P`.
3. It checks that every row is valid for every clustering that separates all far pairs. This check takes one of two forms:
   - **At most 9 vertices:** it enumerates all partitions of the row's vertex set that separate its far pairs.
   - **More than 9 vertices:** the row is accepted only if it is a *star row*. Its `+1` entries are the pairs `{v, t}` for `t` in `T`, and its `-1` entries are a set `R` of pairs inside `T`. The row is accepted if `b <= |T| - |R| - M`:
     - `M = 1` if `R` contains every non-far pair inside `T`;
     - otherwise `M` is the independence number of `(T, R)`.
4. It rounds every `y_i` down to a multiple of `2^-30` and evaluates, in exact integer arithmetic,

   `LB = sum_i b_i y_i + sum_{p in P} ( min(0, c_p - sum_i a_ip y_i) - min(0, c_p) )`,

   where the second sum runs over the pairs that occur in some row.
5. It returns `ceil(LB)`, which is a lower bound on the optimum of the instance: optimal clusterings separate far pairs (Lemma 1), and weak duality applies.
6. Optionally, an archived clustering `LABELS.npz` (key `labels`, one integer per vertex in the numbering above) is evaluated on the same parsed instance.

### Stars as certificates

A pair-disjoint packing of induced stars becomes a valid certificate as follows. For each star with centre `v` and pairwise non-adjacent neighbours `t_1..t_k`, with `k >= 2`, write one row with multiplier 1:
- entries `+1` on the pairs `{v, t_j}`;
- entries `-1` on the pairs `{t_j, t_l}`;
- right-hand side `b = (k-1) - k(k-1)/2`.

Its value is the sum of `k-1` over the stars (Section 3.2 of the paper). `experiments/kapoce/` writes the packings of the KaPoCE branch-and-bound in this way.

## 4. Modes

- `check_certificate.py RAW_FILE CERT.npz`: single file. The format and sign column are taken from the certificate unless the caller passes them. They are printed with the result, which states only that the file, read that way, has the certified bound.
- `check_certificate.py --dir DATA_DIR CERT_DIR OUT.csv` (or `--resume`): every certificate in the directory is checked. `instance_file` must be listed in `DATA_DIR/MANIFEST.sha256` with the same hash. The format and sign column come from `DATA_DIR/INSTANCES.tsv`, never from the certificate.
