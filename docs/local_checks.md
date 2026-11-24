# Local finite checks

These scripts reproduce fixed finite instances used in Chapters 4–6. They do not prove the general theorems or implement a general-purpose solver. Their mathematical inputs, functions, loops, and assertions are preserved from the archived checks.

## Running and dependencies

The repository runner supplies `THESIS_WEB_AUDIT_RESULTS`, the path of an existing output directory. The two scripts write separate JSON files and print a compact result to standard output. Run ordinary Python without `-O`, so the assertions execute.

- `checks/local/local_ch5.py` uses only the Python standard library.
- `checks/local/local_ch346.py` additionally requires NumPy.
- The adapted scripts were tested with Python **3.12.14** and NumPy **2.3.5**.

No thesis source files are needed for the calculations. If the caller explicitly sets `THESIS_WEB_AUDIT_SOURCE`, `local_ch5.py` records SHA-256 identities of `data/chap05.tex`, `data/chap02.tex`, and `data/appendix_w3.tex` relative to that directory and verifies they stay unchanged during execution. Without that optional setting, `source_sha256` is `{}`. The source files are not mathematical inputs. An explicitly provided source directory must contain all three files.

## Chapter 5: six-party, three-level sharing

`local_ch5.py` fixes `q=5` and `L=3` and checks:

| Check | Coverage or expected result |
| --- | --- |
| Authorized-prefix classes | All 64 participant subsets; counts 48, 12, 3, 1 |
| Top-module affine maps | All 63 proper coordinate subsets have full row rank |
| Complete-view representations | 8,000 subset/secret cases; authorized-prefix decoding and hidden-suffix independence |
| Full raw measurement records | For the fixed secret `(1,2,4)`, 390,625 histories map to 15,625 records, each with 25 hidden preimages |
| Fixed-mask outcomes | 125 combinations |
| Resource recovery phase | 3,125 exact modular identities |
| Initial storage and quantum transfers | `N=10`, `C=4`; per-party cut counts agree with the specified schedule |

The view check uses the analytic Bell-module decomposition, the partial-trace argument, affine ranks, and the inverse-Clifford representation. It does not build all density matrices in dimension `5^10`. A full raw record is not necessarily the transcript available to a particular receiver. These finite checks do not establish the general stabilizer-operation lower bounds, arbitrary layer counts, or the general record-refinement lemma.

Output: `local_ch5.json`.

## Chapter 6 phase libraries and Chapter 4 sample counts

`local_ch346.py` fixes `q=5`, `k=3`, `m=1`, and ordered evaluation points `0,1,2,3,4`.

| Encoding phase | LC library size | LU library size |
| --- | ---: | ---: |
| `0` | 25 | 25 |
| `b*c^3` | 125 | 125 |
| `a*c^3+b*c^4` | 625 | 625 |
| `a*c^3+a^4` | 25 | 125 |

LC and LU denote participant-local Clifford and arbitrary unitary operations. Logical phase functions take values in the fifth roots of unity, and an overall phase is removed. The script checks all 25 mask shifts and 625 normalized targets for each encoding: 62,500 pairs, comprising 2,500 feasible and 60,000 excluded pairs.

Each feasible pair is checked by explicit modular back-substitution. Each excluded pair has an integer row satisfying `yB=0` exactly and a nonzero target residual modulo 5. The computed integer kernel rows are not claimed to form a saturated integer basis; every pair is independently resolved by a positive or negative certificate. The output saves the 104 integer kernel rows and **75 negative certificates for three specified targets over all shifts**. Those 75 certificates are not separate saved certificates for all 60,000 excluded pairs.

The script also checks the three original algorithm-A generators on all 125 input rows, giving 375 complete phase identities, including the quadratic remainder of `2*c^3+c^2`.

For the Chapter 4 sample-count table, `epsilon=1/100` and `delta=1/20`. Integer powers check both `N` and `N-1` for gaps `1/10`, `1/4`, `1/2`, and `1/3`; the thresholds are 2,995, 1,197, 598, and 898. This establishes the stated integer thresholds for that bound, not global optimality among all verification strategies.

Output: `local_ch346.json`.

## Nonmathematical portability changes

The public adaptation makes only these changes:

1. In `local_ch5.py`, replace the execution-history docstring with a description of the finite mathematical task.
2. Replace implicit parent-directory discovery and mandatory thesis-source hashing with optional hashing enabled only by `THESIS_WEB_AUDIT_SOURCE`; default to an empty identity map.
3. Remove the obsolete chapter-source line range from the raw-record comment; keep the fixed secret and calculation unchanged.
4. Replace the `identity` result string with a stable mathematical description.
5. Remove two execution-provenance entries from `not_proved`; retain all three mathematical limitations.
6. In `local_ch346.py`, replace the internal task identifier in the sample-count comment with a description of the sample-count table.
7. Replace its `origin` result string with a stable description of the finite checks.

No assertion, mathematical function, fixed parameter, enumeration order, matrix, certificate calculation, or probability comparison is changed. All 24 assertion syntax trees in `local_ch5.py` and all 11 in `local_ch346.py` match the archived versions, as do their mathematical function syntax trees.

A fresh default run matched every mathematical JSON field against the archived results. The only changed top-level fields were `identity`, `not_proved`, and `source_sha256` in `local_ch5.json`, and `origin` in `local_ch346.json`. The `python` field remains runtime metadata and may differ in another environment. Runtime logs and the adaptation audit are kept outside the public repository during preparation; archived reference results are not overwritten.
