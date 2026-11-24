# Finite example checks (E01–E06)

These scripts reproduce finite algebraic and numerical checks used alongside the thesis. They use Python and NumPy, except E01, which uses only the standard library. Each script prints one JSON value to standard output and creates no output files. Run a script from the repository root, for example:

```sh
python checks/examples/E04_ch5_branch.py
```

The repository runner captures this output and validates the reported quantities. Existing assertions are retained, but some scripts report additional quantities without asserting their expected values. In particular, E04 has no internal assertions; a zero exit status alone does not establish that its numerical results passed validation.

## Scope and expected results

| Script | Computation | Expected result |
| --- | --- | --- |
| E01 | Enumerate G, MG, HG and MHG phase functions at five small dimension/system-count pairs; check containments and the G = MG ∩ HG relation; construct null polynomials for composite dimensions up to 60. | State counts listed below; all built-in assertions pass. |
| E02 | Construct every second-layer outcome effect for a three-qutrit example and compare it with the joint spectral effect; compute the equal-weight verification spectrum and sample ceilings. | Spectrum: 0 with multiplicity 16, 0.5 with multiplicity 10, and 1 with multiplicity 1; exact sufficient sample count 598, simplified sufficient count 600. |
| E03 | Compute a rank-two Choi example, a twirl that fixes a selected magic state but erases a selected off-diagonal operator, and classical key/cipher marginals. | Rank 2; nonzero Choi eigenvalues 0.5, 0.5; marginal and twirl errors near zero; all five classical recovery probabilities equal 1. |
| E04 | Evaluate the Choi spectrum and output discrete Wigner function for a two-branch channel on a five-dimensional cubic phase state. | Choi rank 2, trace 1 and spectrum containing 23 zeros and two values 0.5; minimum output Wigner value −√5/25; Wigner sum 1; branch-effect residual near zero. |
| E05 | Check every matrix unit of the logical input against every forbidden subset for two polynomial-code instances, each using one seeded full phase table. | (q,k,n,m) = (5,2,3,1), 100 cross-operator checks; (5,3,5,1), 400 checks; isometry and privacy entrywise errors below 10⁻¹². |
| E06 | Verify finite-field obstruction certificates, polynomial direction ranks, a shift identity and three gate-sequence matrices over F₅. | Direction ranks 2, 3, 4; logical counts 25, 125, 625; all three matrix matches true; 3,125 shift-identity entries vanish. |

E01 expected counts, in the order G, MG, HG, MHG:

| d | n | Counts |
| --- | --- | --- |
| 2 | 2 | 2, 2, 8, 8 |
| 3 | 2 | 3, 81, 27, 6,561 |
| 4 | 1 | 1, 1, 4, 16 |
| 4 | 2 | 4, 64, 64, 16,384 |
| 6 | 1 | 1, 1, 6, 18 |

The composite dimensions checked by E01 are 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 34, 35, 36, 38, 39, 40, 42, 44, 45, 46, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, and 60.

## Quantities requiring output validation

The following reported values are not all covered by the scripts' own assertions. The runner should compare them with the expected values, using exact equality for integers and booleans and an absolute tolerance of 10⁻¹² for zero residuals. A relative tolerance of 10⁻¹⁰ and absolute tolerance of 10⁻¹² accommodates floating-point library variation in the other reported numbers. E02 serializes spectrum values as dictionary keys; normalize `-0.0` and `0.0` to the same numeric key before comparison.

- **E02:** `allpass_bound_N598 = 0.04991169123005856`, `allpass_bound_N597 = 0.05016250374880256`, and `correlated_counterexample_allpass_lower_bound = 0.99`. The last field is the scalar `1 − epsilon`; the script does not construct or simulate a correlated source. The internal assertions already check the threshold crossing at 598.
- **E03:** `nonzero_choi_eigenvalues = [0.5, 0.5]`; `choi_reference_marginal_error`, `unrecorded_twirl_magic_fixed_error`, `key_view_max_difference` and `cipher_view_max_difference` must be near zero. The Choi rank, off-diagonal twirl residual and recovery probabilities have internal assertions as well.
- **E04:** validate every computed diagnostic: `choi_eigenvalues = [0] * 23 + [0.5, 0.5]`, `choi_rank = 2`, `choi_trace = 1`, `output_wigner_min = -0.08944271909999159`, `output_wigner_sum = 1`, and `max_branch_effect_error ≈ 0`. The `q = 5` and `ancilla_probabilities = [0.5, 0.5, 0, 0, 0]` fields declare the chosen example; the probabilities are written directly in the result dictionary rather than independently computed.
- **E05:** the two `cross_operator_checks` values must be 100 and 400. Seeds 20261027 and 20261028 identify the selected phase tables; only those tables are tested.
- **E06:** `rows = 125`, `B_shape = [125, 25]`, `B_nonzeros = 625`, `D27_yB_integer_max = 30`, `polynomial_direction_ranks = [2, 3, 4]`, `logical_counts_from_direction_ranks = [25, 125, 625]`, and `light_shift_identity_checks = 3125`. `D23_cyclic_checks = 25` labels the number of cases in the preceding asserted loop. Both modular products `D27_yB_zero_mod5` and `D27_yJ_zero_mod5` are true, `D27_y_minus_delta_h_mod5 = 1`, while `D27_yB_exact_integer_zero` is false. Thus a modular annihilation certificate is distinguished from an integer one. The `D10`, `D11`, `D12`, `D23` and `D27` prefixes are stable output identifiers and are not promises of equation numbering in another thesis edition.

E01's enumerated counts also require comparison with the count table above; its built-in assertions check divisibility and set relations, rather than those particular cardinalities.

## Optional source integrity metadata

E02 performs its mathematical computation without thesis source files. By default, `source_sha256` is `{}`. If the caller explicitly sets `THESIS_AUDIT_SOURCE_ROOT`, E02 hashes the six relative files listed in its `names` variable before and after computation and requires identical hashes. A missing file then raises an error. `source_hashes_unchanged` records this equality; with an empty hash dictionary, no thesis source was examined. Hashes are provenance metadata and are excluded from numerical reference comparisons.

## Interpretation

These are reproducible checks for the stated finite examples. E01 does not enumerate all dimensions; E05 does not exhaust all phase functions or multisystem logical inputs. E03 and E04 isolate specific channel and secrecy distinctions. The scripts retain the thesis's general mathematical proofs as separate arguments. NumPy roundoff can change residuals near machine precision without changing the checked conclusion.
