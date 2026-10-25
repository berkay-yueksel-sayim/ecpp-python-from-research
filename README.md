# ecpp-python-from-research

**Prove that a large number is prime, and check the proof: elliptic curve primality proving (ECPP) in pure Python.**

ECPP turns "this number passed every probabilistic test" into a proof: a certificate, a short chain of elliptic curves and
points, that anyone can check independently, without trusting the program that found it. This package contains

- **a certificate checker** that needs nothing but the Python standard library and reads the certificates that PARI/GP's
  `primecert` produces with its default flag, and
- **a prover** that builds such certificates (Atkin–Morain method with complex multiplication; needs `mpmath`, uses
  `gmpy2` if installed).

The code was written for a paper series on consecutive powerful numbers, where it checked every elliptic curve
certificate the papers rely on. The research version and its certificates are kept in [`research/`](research/), as
byte-identical copies of the files in the data repository of the series (with a README added). Before this release, a
review found a gap in the curve arithmetic of the research checker; this package corrects it (see
[A correction](#a-correction)).

## Quick start

Requires Python 3.10 or later.

```bash
git clone https://github.com/berkay-yueksel-sayim/ecpp-python-from-research
cd ecpp-python-from-research
pip install ".[prove]"            # ".[prove,fast]" also installs gmpy2 for a faster prover
ecpp-research verify research/ergebnisse/w164_zertifikate/cert_6.txt
ecpp-research prove 1000000000000000000000000000000000000000000000000000000000007 --out cert.txt
ecpp-research verify cert.txt --number 1000000000000000000000000000000000000000000000000000000000007
```

The distribution is named `ecpp-research`, the module `ecpp_research`; `python -m ecpp_research …` works like
`ecpp-research …`. `verify` prints the number the certificate proves (shortened if long; `--number` checks that it is
the one you expect) and exits with 0 if the certificate is valid, 1 if it is not, and 2 if the file cannot be read or
is not a certificate in the expected format. `prove` (options `--out`, `--time-limit` in seconds, default 3600, `--seed`,
`-v`) exits with 0 when it found and verified a certificate, 1 when it did not (N at most 2^64, not a probable prime, or no
certificate within the limits), and 2 when N is not an integer.

```python
from ecpp_research import check_certificate, load
result = check_certificate(load("cert.txt"))
print(result.valid, result.steps, result.reason)

from ecpp_research.prove import prove
steps = prove(2**127 - 1)            # a list of steps (N, t, s, a, (x, y)), as tuples of ints
```

`check_certificate` takes a list of such steps; for anything else (a step of the wrong length, a missing point, a value
that is not an integer) it returns an invalid result with the reason instead of raising.

## Where it comes from, and what it supports

The paper series *Consecutive Powerful Numbers and Pell Equations: Triples, Pairs, and Further Observations* (B. Y. Sayim,
2026) consists of a volume and three parts:

| | Zenodo (all versions) |
|---|---|
| Volume | [10.5281/zenodo.23127546](https://doi.org/10.5281/zenodo.23127546) |
| Part I — *Triples of Consecutive Powerful Numbers: Pell Addresses, a Height Inequality, and Witnesses* | [10.5281/zenodo.23127702](https://doi.org/10.5281/zenodo.23127702) |
| Part II — *Pairs of Consecutive and Near-Consecutive Powerful Numbers: Pell Families, Towers, and Localization* | [10.5281/zenodo.23127704](https://doi.org/10.5281/zenodo.23127704) |
| Part III — *Further Observations on Powerful Numbers: Gaps, Cyclotomic Values, and Independence Tests* | [10.5281/zenodo.23127706](https://doi.org/10.5281/zenodo.23127706) |

Data, scripts and result files of the series are in the data repository
[erdos-364-365-powerful-numbers-pell](https://github.com/berkay-yueksel-sayim/erdos-364-365-powerful-numbers-pell).

A number is powerful if every prime that divides it divides it at least twice, so one prime that divides a number exactly
once shows that the number is not powerful. The papers use such primes as witnesses, together with complete
factorizations, and prove the primes they use; 25 of them are proved prime by elliptic curve certificates (page numbers
of v1.0; the PDFs of v1.0 and v1.0.1 are identical):

- **Part I, Proposition I.6.8** (stated on p. 26; the certificates are described in its proof on p. 28): the 6 largest
  prime factors used in the complete factorizations and in the algebraic split of building blocks, with 111 to 477 digits.
- **Part II, Proposition II.7.2 (b)** (p. 21): 19 primes with 48 to 1687 digits, the prime cofactors in (b) that lie above
  `psi_12` (the other 5 lie below it and are proved by the deterministic test).

PARI/GP 2.17.2 generated these 25 certificates, and the research version of this checker verified each of them step by
step. Our own prover, written as a second route, found certificates for 23 of the 25 primes (all 6 of Part I, 17 of the
19 of Part II, up to 693 digits), and PARI/GP's verifier `primecertisvalid` accepted all 23. The two largest primes of
Part II (1191 and 1687 digits) have only the PARI/GP certificates.

All 48 certificates for these primes are accepted by the corrected checker of this package and by PARI/GP's
`primecertisvalid` ([`checks/pari_all_25_output.txt`](checks/pari_all_25_output.txt) for the 25 PARI/GP certificates).

This repository makes that code a standalone package: English names, a command line, tests, and this documentation.
Part III uses no elliptic curve certificates; the volume collects the three parts. The certificate files are numbered
by computation job (`cert_<i>.txt`, with the numbers listed in `w164_auftrag.json` and `w165_auftrag.json`); the papers
do not number these primes.

## A correction

The research checker (`research/skripte/w164_ecpp_pruefer_2026-10-02.py`) and the first draft of this package decided
the case distinctions of the curve arithmetic (equal points, opposite points) modulo N only, and treated degenerate
points as the point at infinity. For a composite N a case can hold modulo one prime divisor and fail modulo another; the
arithmetic then gave wrong results modulo that divisor without noticing. A review before this release (23 October 2026)
built certificates for composite numbers that the research checker accepts, for example one "proving"
10196873 = 1637 · 6229 prime ([`tests/test_soundness.py`](tests/test_soundness.py)).

The checker of this package now uses an arithmetic that refuses every such case: an element is used only when it is 0
modulo N or invertible modulo N, so every computed point is the true point modulo every prime divisor of N, or the step is
rejected (argument in [`ecpp_research/curves.py`](ecpp_research/curves.py)). It accepts all 50 certificates in
`research/` and rejects the false ones ([`tools/test_run_output.txt`](tools/test_run_output.txt),
[`tools/benchmark_output.txt`](tools/benchmark_output.txt)).

The results of the papers are not affected: their certificates are valid, as the corrected checker shows. What changes is
the statement that the research program alone verified them; that program could also have accepted invalid certificates
of a special form.

In the same review, a comparison with PARI/GP found that PARI/GP's `primecertisvalid` (version 2.17.2) also accepts some
certificates for composite numbers: 6 of 1203 one-step certificates for products of two primes
([`checks/compare_with_pari_output.txt`](checks/compare_with_pari_output.txt), with the exact recomputation of each).
Version 2.19.0 (Debian package `pari-gp` 2.19.0-1), the newest at the time of writing, gives the same 2403 verdicts, byte
for byte (`checks/compare_pari.txt`). [`checks/pari_minimal_reproduction.gp`](checks/pari_minimal_reproduction.gp)
reproduces the six cases with PARI/GP alone, without this package: for each it prints PARI/GP's verdict, the factorization
of N, and q*(s*P) modulo each prime factor, which is not the point at infinity modulo any of them (output for 2.17.2 and
2.19.0 in [`checks/pari_minimal_reproduction_output.txt`](checks/pari_minimal_reproduction_output.txt)). We reported this
to the PARI/GP developers on 25 October 2026 (bug #2723 in the PARI/GP bug tracker). PARI/GP's agreement is therefore a second opinion, not a proof.

## Scope and evidence

*Stage: algebraic = a statement about all numbers; finite-size numeric = computed for the numbers named, in exact integer
arithmetic unless stated. Label: what kind of statement it is.*

| Statement | Stage | Label |
|---|---|---|
| If a certificate satisfies conditions 1–5 below modulo every prime divisor of N, its first number N is prime. | algebraic | theorem (Atkin–Morain 1993, Theorem 5.2 and Corollary 5.1, after Goldwasser–Kilian 1986; for the last link, Sorenson–Webster 2017, Theorem 1.1) |
| Every point computed by the checker's arithmetic is the true point modulo every prime divisor of N, or the step is rejected. | algebraic | theorem for the arithmetic as written (argument in `ecpp_research/curves.py`); the code is tested, not formally verified |
| On 4000 random cases with N = p·r, every result of the arithmetic equals the exact result modulo p and modulo r, or is refused (1035 compared, 0 wrong); the arithmetic of the research checker gives a wrong result modulo a prime divisor in 2122 of the same 4000 cases (2155 of the 8000 pairs of case and prime divisor, which is what the test counts). | finite-size numeric | exact finite |
| The checker accepts all 50 certificates in `research/` (the 25 PARI/GP certificates the papers rely on, the 23 of our own prover, two controls for nextprime(10^40)) and rejects the three false certificates for composite numbers that the research checker accepts. | finite-size numeric | exact finite |
| On 2403 one-step certificates (240 valid and 960 changed ones for primes with 4 or 5 digits, 1203 for composite numbers), this checker and PARI/GP's `primecertisvalid` give the same verdict, except for 6 certificates for composite numbers that PARI/GP accepts and this checker rejects. | finite-size numeric | exact finite (`checks/compare_with_pari_output.txt`) |
| The 25 primes of Parts I and II (48 to 1687 digits) are prime. | finite-size numeric | computer-assisted theorem: exact integer arithmetic, each certificate accepted by the corrected checker (rows above); PARI/GP's verifier accepts them as well |
| In every test and benchmark run of this package the prover found a certificate: 9 random primes of 30, 45 and 60 digits and the next primes after 10^30, 10^35 and 10^40 (tests), the next primes after 10^40, 10^60, 10^100, 10^150 and 10^200 (benchmark), and the six numbers of Part I with 111 to 477 digits (slow test). | finite-size numeric | exact finite for these numbers; running times are single runs on one machine (table below) |
| The Hilbert class polynomials are computed in multiprecision floating point (`mpmath`), on two routes, and rounded. | finite-size numeric | numerical evidence; no proof is needed, because a step built from them is kept only if the checker accepts it, so a wrong polynomial can cost time but cannot produce a false proof (Atkin–Morain 1993, p. 41) |

## How a certificate is checked

A certificate is a list of steps `[N, t, s, a, [x, y]]`, the format of PARI/GP's `primecert`. For each step, with
`m = N + 1 - t`, `q = m / s` and `b = y^2 - x^3 - a*x (mod N)`, the checker requires

1. `gcd(N, 6) = 1` (so every prime divisor is at least 5 and the short Weierstrass form applies) and `t^2 < 4N`;
2. `s` divides `m` and `q > (floor(N^(1/4)) + 2)^2`, which implies `q > (N^(1/4) + 1)^2`;
3. the curve `y^2 = x^3 + a*x + b` is nonsingular modulo every prime divisor of `N`: `gcd(4a^3 + 27b^2, N) = 1`;
4. `Q = s*P` is not the point at infinity modulo any prime divisor of `N`, and `q*Q` is, modulo every prime divisor. Both
   are computed in affine coordinates modulo `N`, where each case of the group law (equal points, opposite points, a
   vanishing denominator) is decided only when it holds modulo every prime divisor at once: a quantity is used only if it
   is 0 modulo `N` or invertible modulo `N`, and anything else rejects the step;
5. the next step starts with `N = q`, and the last `q` is at most `2^64` and passes the Miller–Rabin test with the first
   twelve prime bases, which is deterministic below `psi_12 = 318665857834031151167461` (Sorenson–Webster 2017,
   Theorem 1.1).

Conditions 1 to 4 contain the hypotheses of Atkin–Morain 1993, Theorem 5.2 and Corollary 5.1 (pp. 37–38), which they
attribute to Goldwasser and Kilian; the test `t^2 < 4N` is an additional consistency check (it holds for every prime N by
Hasse's bound) and is not needed for the conclusion. If all conditions hold, every prime divisor of `N` exceeds
`sqrt(N)`, so `N` is prime once `q` is. Certificates are parsed with `ast.literal_eval`, which evaluates no code.

## How the prover works

For a probable prime `N > 2^64` the prover looks for a discriminant `-D` with `4N = A^2 + D*B^2` (Cornacchia), takes the
candidates `m = N + 1 ± A`, removes the primes below `3 * 10^6` from `m`, and keeps those where the rest `q` is a probable
prime above the bound (Atkin–Morain 1993, Section 5.3). From a root of the Hilbert class polynomial `H_D` modulo `N` it
builds a curve with that `j`-invariant (Section 8.6.2, Equation (39)) and, by a scaling that needs no square roots
(Section 8.6.3), a point on that curve or on its quadratic twist, which has the same `j`-invariant; it keeps the step only
if the checker accepts it. Then it continues with `q`. Discriminants come in three lists of growing size (class number up
to 8, 14 and 24); when no step is found, the prover goes back one step and excludes that `q`, at most 40 times.

The prover only searches. Every step it returns has passed the checker, and the finished chain is checked again, so an
error in the prover can cost time but cannot produce a false proof.

## Running times

Single runs of [`tools/benchmark.py`](tools/benchmark.py) on one laptop (Intel Core processor; the code uses one core),
Python 3.12 on Windows 11, `gmpy2` installed for the prover, on mains power, while one other computation ran on the same
machine ([`tools/benchmark_output.txt`](tools/benchmark_output.txt), which lists all 50 certificates; run with version
0.1.0, whose arithmetic version 0.1.1 keeps unchanged). Runs of the same
checker on battery power took up to about half as long again. Take the numbers as an order of magnitude.

| Task | Digits | Steps | Time |
|---|---|---|---|
| check the PARI/GP certificate of the largest prime of Part I | 477 | 67 | 15 s |
| check a PARI/GP certificate of Part II | 693 | 93 | 50 s |
| check a PARI/GP certificate of Part II | 1191 | 157 | 7 min (419 s) |
| check the PARI/GP certificate of the largest prime of Part II | 1687 | 184 | 22 min (1311 s) |
| prove the next prime after 10^100 | 101 | 7 | 1.6 s |
| prove the next prime after 10^200 | 201 | 16 | 7.1 s |

The checker needs one modular inversion for every group operation; that is the price of deciding every case modulo all
prime divisors at once (see [A correction](#a-correction)). The first call of `prove` in a session is slower (1.9 s for
the 41-digit prime, then 0.3 s for the 61-digit one): it computes the product of all primes below `3 * 10^6` once and
caches the class polynomials it computes.

## Limits

- Pure Python. Checking the 1687-digit certificate took 22 minutes (see Running times). The prover is meant for numbers up to several hundred digits. In the research runs, PARI/GP
  2.17.2 generated the certificate of the 477-digit prime in 8 seconds, where our research prover needed 11 minutes (673 s,
  in a second run with one backtrack, after a first run found no step at 388 digits). PARI/GP generated the certificates
  of the primes of 1191 and 1687 digits in about 4.5 and 13 minutes; our research prover produced no certificate for them
  (for 1687 digits it found no first step within 344 seconds, and the run for 1191 digits left no result). These PARI/GP
  times show the order of magnitude, not a controlled comparison (`research/ergebnisse/w164_zertifikate/zeit.txt` and
  `w165_zertifikate/zeit.txt`, wall-clock milliseconds; `research/ecpp/ergebnisse/e5_teil1_log.txt` and `e5_teil2_log.txt`).
- PARI/GP's `primecertisvalid` is much faster (77 seconds for the 1687-digit certificate, `checks/pari_all_25_output.txt`)
  but, in our comparison, accepted some certificates for composite numbers (see A correction).
- The implementation is tested, not formally verified.
- It is a mathematical tool, not written for cryptographic products: the arithmetic is not constant-time and has no
  protection against side channels. Patents exist on validating or generating primes for use in cryptographic protocols
  (for example US 11,251,953 B2 and US 11,411,732 B2); whoever builds this code into such a product should check them for
  that use.
- Certificates are parsed with `ast.literal_eval`, which evaluates no code but is not hardened against very large inputs;
  limit the file size when checking files from untrusted sources. On Python 3.11 and later, integers above 4300 digits are
  read and written with Python's integer-string limit lifted while the certificate is parsed or written, and for the whole
  run of the command line.
- The prover uses discriminants with class number at most 24 and `D < 200000`.
- Numbers up to `2^64` need no certificate: use `is_prime_small`, which is deterministic there.

## Tests

```bash
pip install ".[test]"
pytest                      # about two minutes
ECPP_SLOW=1 pytest          # adds the 8 certificates above 500 digits, the class polynomials for D < 2000, and the re-proof of Part I (about an hour on the benchmark machine, tools/test_run_output.txt)
```

The tests carry over the controls of the research scripts (Jacobi symbol against Euler's criterion, square roots,
Cornacchia against exhaustive search, the class numbers of Atkin–Morain's Table 1, the class polynomials of class number
1 and `H_23`, the reference values of their pp. 59–60), check every certificate in `research/` with both the package
checker and the research checker, change the first step of the small certificates at random (800 seeded changes), and
run the prover on random primes of 30, 45 and 60 digits. [`tests/test_soundness.py`](tests/test_soundness.py) checks the
arithmetic against exact computation modulo each prime divisor (with the research arithmetic as positive control) and
rejects the three false certificates; malformed input (13 kinds) must give an invalid result, not an exception.
[`checks/`](checks/) holds the comparison with PARI/GP and the minimal reproduction, which need PARI/GP and are not part
of `pytest`. The automatic tests on GitHub run on Python 3.10 to 3.13, and once with `gmpy2` installed.

## Related software

Established ECPP programs include PARI/GP (`primecert`, `primecertisvalid`), which generated the certificates of the
papers, and Primo; implementations exist in several other languages. This package was written from the published papers
listed below and, for the certificate format, from the PARI/GP documentation of `primecert` and `primecertisvalid`,
without reading the source code of any other implementation; outputs were compared, never code.

## References

- S. Goldwasser and J. Kilian, *Almost all primes can be quickly certified*, Proceedings of the Eighteenth Annual ACM
  Symposium on Theory of Computing (STOC '86), 316–329, 1986. [doi:10.1145/12130.12162](https://doi.org/10.1145/12130.12162)
- A. O. L. Atkin and F. Morain, *Elliptic curves and primality proving*, Math. Comp. 61 (1993), no. 203, 29–68.
  [doi:10.1090/S0025-5718-1993-1199989-X](https://doi.org/10.1090/S0025-5718-1993-1199989-X)
- J. Sorenson and J. Webster, *Strong pseudoprimes to twelve prime bases*, Math. Comp. 86 (2017), no. 304, 985–1003.
  [doi:10.1090/mcom/3134](https://doi.org/10.1090/mcom/3134); preprint [arXiv:1509.00864](https://arxiv.org/abs/1509.00864),
  where the bound `psi_12` is Theorem 1.1.
- The PARI Group, *PARI/GP*, version 2.17.2, Univ. Bordeaux; documentation of the functions `primecert` and
  `primecertisvalid`.

## License and citation

The code and the documentation are licensed under the Apache License, Version 2.0 (`LICENSE`, `NOTICE`). The certificates,
result files and PARI/GP job files in `research/ergebnisse/` and `research/ecpp/ergebnisse/`, and the run outputs in
`tools/` and `checks/`, are licensed under CC BY 4.0 (`LICENSE-data`), as in the data repository.

If you use the checker, the prover or the certificates, please cite the software,
[doi:10.5281/zenodo.23585827](https://doi.org/10.5281/zenodo.23585827) (all versions; each version also has its own DOI
on Zenodo), and the volume, [doi:10.5281/zenodo.23127546](https://doi.org/10.5281/zenodo.23127546).

Generative AI tools were used in writing this code and its documentation, as in the papers of the series, where the use is
disclosed.
