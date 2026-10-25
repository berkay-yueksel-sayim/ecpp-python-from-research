# Changelog

## 0.1.1

- `check_certificate` and `check_step` return an invalid result with a reason for input that is not a list of steps
  `(N, t, s, a, (x, y))` with integers (a step of the wrong length, a missing point, a float or a string), instead of
  raising an exception. No valid or invalid certificate is judged differently; the curve arithmetic is unchanged.
  New tests for 13 kinds of malformed input.
- `checks/pari_minimal_reproduction.gp` and its output: the six certificates for composite numbers that PARI/GP's
  `primecertisvalid` accepts, reproduced with PARI/GP alone (2.17.2 and 2.19.0 give the same output). The finding is
  bug #2723 in the PARI/GP bug tracker.
- The automatic tests run on Python 3.10, 3.11, 3.12 and 3.13, and on 3.12 also with `gmpy2`.
- The README asks to cite the software by its concept DOI (all versions). `tools/benchmark_output.txt` is the run of
  version 0.1.0; `tools/test_run_output.txt` is a new run of the complete test suite.

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
