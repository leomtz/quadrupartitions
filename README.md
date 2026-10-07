# Readme

Simulated annealing and exact verification for orthogonal partition counterexamples.

A small Python/Jupyter companion to [**Counterexamples and symmetry for uneven
orthogonal mass partitions in the plane**](https://arxiv.org/abs/2609.16757),
by Leonardo Martínez-Sandoval (arXiv:2609.16757, 2026).

The notebook explains a polar simulated-annealing search, generates a fresh
centrally symmetric configuration of 96 integer points, and certifies that it
has no orthogonal partition with cyclic counts `(8, 8, 40, 40)`, even allowing
points on the cutting lines to be assigned to adjacent sectors. It separately
checks the configuration printed in the article.

## Run the notebook

Use Python **3.11 or newer** (tested on 3.12.14). From the repository directory:

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m jupyterlab
```

Open `annealing.ipynb`, select the environment's Python kernel, and choose
**Restart Kernel → Run All**. Keep the notebook beside the Python modules.
All notebook explanations are in English, with rendered mathematical formulas.

The search starts from a fresh irregular ring and uses the calibrated seed `5`.
In the recorded environment it reached zero energy after **16,947 proposals**,
in about **4 seconds**; both exact checkers together took about **2 seconds**
per example. Compilation, figures, and notebook startup add overhead.
The default search budget is **240 seconds**, and it stops early on success.
Other seeds or environments can have different outcomes and runtimes. If the
budget expires, the notebook raises an explicit error; it never silently loads
a saved counterexample in place of a failed search.

To execute headlessly (allowing time for compilation, search, and certificates):

```bash
python -m nbconvert --to notebook --execute annealing.ipynb \
  --ExecutePreprocessor.timeout=900 --output annealing.executed.ipynb
```

## Verify the article example without scientific Python dependencies

The exact verifiers use only the Python standard library:

```bash
python exact_verify.py data/paper_example.json --full
python slope_verify.py data/paper_example.json
# Run both, plus exact geometric checks:
python verification.py data/paper_example.json
```

Expected: `counterexample_certified`, 9,216 open chambers, overlap-defect
histogram `{1: 8814, 2: 402}`, and zero weak solution events. The first point is
`[6501822, -2132732]`. The checked-in `data/python_example.json` is the separate
calibrated Python result; it can be verified in the same way.

## Mathematical and computational scheme

Write `P = {p_1, ..., p_48, -p_1, ..., -p_48}`. For a cap normal `n`, choose
the median half-set `A = Top_48(Jn)` and cap `B = Top_16(n)`. The target
partition occurs exactly when `q = |A intersect B| = 8`.

The search minimizes the squared shortfall below `q = 9` over its floating-point
directional event sweep. It moves one representative in polar coordinates,
accepts improvements and occasional worse moves using a Metropolis rule,
and decreases temperature and move size over a cooling cycle. Restarts begin
with new irregular rings. NumPy and Numba accelerate Python source; no C++
code is needed to run this repository.

Search and certification are separate. Search energy zero alone is not a proof.
After rounding representatives to integers and adding their negatives:

- `exact_verify.py` checks all integer event rays and intervening chambers,
  and the overlap ranges attainable at tied cutoffs.
- `slope_verify.py` independently checks rational-slope charts and uses a
  dynamic program for critical-direction assignments.
- `verification.py` checks agreement, central symmetry, distinctness, and
  all 142,880 triple determinants.

## How to cite

If you use the code or configurations in this repository, please cite
**quadrupartitions**. If you use the mathematical results or construction,
please also cite the accompanying article.

### Repository

Leonardo Martínez-Sandoval. *quadrupartitions: Simulated annealing and exact
verification for orthogonal partitions*. GitHub repository, 2026.
[https://github.com/leomtz/quadrupartitions](https://github.com/leomtz/quadrupartitions).

```bibtex
@misc{martinezsandoval2026quadrupartitions,
  author       = {Mart{\'i}nez-Sandoval, Leonardo},
  title        = {{quadrupartitions}: Simulated annealing and exact verification
                  for orthogonal partitions},
  year         = {2026},
  howpublished = {GitHub repository},
  url          = {https://github.com/leomtz/quadrupartitions}
}
```

### Article

Leonardo Martínez-Sandoval. *Counterexamples and symmetry for uneven orthogonal
mass partitions in the plane*. arXiv:2609.16757, 2026.
[https://arxiv.org/abs/2609.16757](https://arxiv.org/abs/2609.16757).

```bibtex
@misc{martinezsandoval2026orthogonalpartitions,
  author        = {Mart{\'i}nez-Sandoval, Leonardo},
  title         = {Counterexamples and symmetry for uneven orthogonal mass
                   partitions in the plane},
  year          = {2026},
  eprint        = {2609.16757},
  archivePrefix = {arXiv},
  url           = {https://arxiv.org/abs/2609.16757}
}
```

## Files and outputs

| File | Purpose |
|---|---|
| `annealing.ipynb` | Executed tutorial with formulas, code, plots, and certificates |
| `annealing.py` | Fresh-ring initialization, event sweep, polar moves, and annealing |
| `exact_verify.py`, `slope_verify.py` | Independent exact verifiers, preserved from the historical audit |
| `verification.py` | Readable combined report and geometric checks |
| `data/paper_example.json` | Exact article coordinates |
| `data/python_example.json` | Exact coordinates of the calibrated Python result |
| `PROVENANCE.md` | Historical C++ run and Python calibration record |
| `outputs/` | Generated coordinates, reports, figures, and environment metadata |

Saved notebook outputs document one real execution. Running it again performs
a new search and overwrites the generated output files. The new candidate's
`status` stays `UNCERTIFIED` in its search JSON; the separate certificate is
authoritative. Coordinate hashes refer to the ordered list, not file bytes.

For quick checks of the sweep against independent full projection orders,
boundary assignments, and a positive partition example:

```bash
python -m unittest -v test_verification
```
