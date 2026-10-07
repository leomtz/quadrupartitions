# Provenance and reproducibility

## Configuration in the article

`data/paper_example.json` is the historical
`candidate_default_seed304_N96_k8.json` from the package
`reproduction_N96_k8_2026-09-14.zip`. All 96 points were compared against the
table and central-symmetry rule in the supplied `arxiv-V1.tex`; they agree.
Its first representative is `(6501822, -2132732)`.

The experiment dated 2026-09-14 compiled the user-supplied `anneal.cpp`
unchanged and ran:

```bash
g++ -O3 -std=c++17 anneal.cpp -o anneal
mkdir -p discrete_small_k
./anneal 48 8 20000 304 sweep
```

It reached zero energy at restart index 5, step index 15281. The experiment
was described as a new reproduction, distinct from an earlier configuration;
the newly reproduced configuration is the one in the present article.

The exact audit reported 9,216 chambers, defect histogram
`{1: 8814, 2: 402}`, no weak solution events, and no collinear triples.
Ordered-coordinate SHA-256:

```text
2b314645c5a8fd3f44b4334f539db860d658a5476ab2dff662e800ceabb18ed9
```

The two standalone verifiers here preserve the historical `exact_verify.py`
and `slope_verify.py` source. They are independent of the floating search.

## Python search calibration (2026-10-07)

The Python version preserves the search mechanism, not the C++ random stream.
Numba compiles Python/NumPy functions; it does not invoke the historical C++
program. Initial representatives are fresh irregular rings. No counterexample
file is read by `search()`.

The following seeds were explored during calibration:

| Seed | Restart budget | Proposals | Best energy |
|---:|---:|---:|---:|
| 304 | 8 | 160,000 | 11 |
| 1 | 6 | 120,000 | 117 |
| 2 | 6 | 120,000 | 23 |
| 3 | 6 | 120,000 | 57 |
| 4 | 6 | 120,000 | 30 |
| 5 | 6 | 16,947 | 0 |

Each restart allowed 20,000 proposals, with progress chunks of 1,000.
Seed 5 succeeded on restart index 0 in about 3.5 seconds of search.
This is a selected demonstration seed, not an unbiased success-rate study.
The rounded result was certified by both verifiers: 9,216 chambers,
histogram `{1: 8948, 2: 268}`, zero weak solution events, 96 distinct points,
central symmetry, and zero collinear triples among 142,880 triples.
Its first representative is `(9535992, -2805451)`; its ordered-coordinate hash is:

```text
7745d0f3e70dfc87d6bf0b35ff88af1dc9dc4d667105b9de20e6ae209f57d207
```

The checked-in `data/python_example.json` records this calibrated result.
The notebook generates it again rather than loading it as a starting state.

Calibration used Python 3.12.14, NumPy 2.3.5, Numba 0.68.0, and Matplotlib 3.10.8.
The notebook saves dependency versions, parameters, attempts, coordinate hashes,
and source hashes in `outputs/run_metadata.json`.

The Python routine reseeds its compiled generator at chunk boundaries using
deterministic seeds supplied by NumPy's `Generator`. Consequently chunk size
is part of the experiment. Compiler, platform, dependency versions, or a
wall-clock cutoff can affect the trajectory. Exact verification of a fixed
integer file does not depend on these search details.

## Limits of the computation

The floating search can misorder near-coincident events, and rounding can
change the geometry. Every assertion of a counterexample is therefore made
only after exact verification of the integer output. Failure to find one
within a finite budget proves nothing. Neither minimality of N=96 nor
minimality of k=8 is asserted.
