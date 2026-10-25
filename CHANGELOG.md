# Changelog

## 0.1.0

First public version.

- `ecpp_research.verify`: certificate checker for the PARI/GP `primecert` format, standard library only; the conditions
  of the research checker that verified the certificates of the paper series, in the same order.
- `ecpp_research.curves`: affine group arithmetic modulo N that decides each case (point at infinity, equal points)
  only when the decision holds modulo every prime divisor of N, and otherwise refuses the step. The research checker
  decides these cases modulo N only and accepts certificates for some composite numbers (`research/README.md`); the
  certificates of the papers are not affected. `tests/test_soundness.py` holds three such certificates and a property
  test of the arithmetic against each prime divisor.
- `ecpp_research.prove`: Atkin–Morain prover (Cornacchia, Hilbert class polynomials computed on two routes, best of
  several candidates, three discriminant levels, backtracking one step at a time, up to 40 times); every step is
  checked before it is kept.
- Command line `ecpp-research verify` / `ecpp-research prove` (also `python -m ecpp_research`); exit code 2 for any
  text that is not a certificate in the expected format.
- Certificates of numbers above 4300 digits are read and written with Python's integer-string limit lifted for the
  duration of the call. Requires Python 3.10 or later.
- Tests carried over from the controls of the research scripts, cross-checks against the unchanged research checker,
  and `tools/benchmark.py` with its output.
- `research/`: the research version and all its certificates, unchanged.
- `checks/pari_all_25.gp` and its output: PARI/GP's `primecertisvalid` accepts all 25 PARI/GP certificates of the papers.
- `checks/compare_with_pari.py` and `.gp` with their outputs: the checker and PARI/GP 2.17.2's `primecertisvalid` on
  2403 certificates (valid ones, corrupted ones, and one-step certificates for composite numbers).
