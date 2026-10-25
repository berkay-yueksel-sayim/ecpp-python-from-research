# SPDX-License-Identifier: Apache-2.0
"""Soundness: the checker must never accept a certificate for a composite number.

Found in review on 23 October 2026: the curve arithmetic of the research checker (and of this package before the fix)
decided case distinctions modulo N only and treated degenerate points as the point at infinity, so it accepted the
three certificates below for composite numbers. The second test checks the property that rules this out: every result
of ``curves`` is the true result modulo every prime divisor of N, or the operation refuses (``curves.Degenerate``).
The same harness applied to the old arithmetic finds violations; that is its positive control."""
import math
import random

from ecpp_research import check_certificate, curves

COMPOSITE_CERTIFICATES = [
    [(10196873, 71, 1, 1021388, (4608684, 5613837))],
    [(10210063, 3265, 1, 586389, (4233552, 698090))],
    [(10270709, -1301, 1, 5725870, (4343939, 1328555))],
]


def smallest_factor(n):
    return next((d for d in range(2, math.isqrt(n) + 1) if n % d == 0), n)


def test_known_false_certificates_are_rejected(w164):
    for cert in COMPOSITE_CERTIFICATES:
        N = cert[0][0]
        assert smallest_factor(N) < N                                  # N is composite
        assert w164['pruefe_zertifikat']([(*s[:4], list(s[4])) for s in cert]) is True  # the research checker's gap
        result = check_certificate(cert)
        assert not result.valid and 'degenerate' in result.reason


def _ref_add(P, Q, a, p):
    """The group law over F_p, with every case."""
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def _ref_multiply(k, P, a, p):
    R = None
    for bit in bin(k)[2:]:
        R = _ref_add(R, R, a, p)
        if bit == '1':
            R = _ref_add(R, P, a, p)
    return R


def _cases(seed, count):
    """Random (N = p*r, a, P, k) with small primes p, r, so that cases holding modulo p only are frequent."""
    rng = random.Random(seed)
    primes = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    out = []
    while len(out) < count:
        p, r = rng.sample(primes, 2)
        N = p * r
        a, x, y = rng.randrange(N), rng.randrange(N), rng.randrange(N)
        b = (y * y - x * x * x - a * x) % N
        if any((4 * a ** 3 + 27 * b * b) % f == 0 for f in (p, r)):
            continue
        out.append((p, r, N, a, (x, y), rng.randrange(1, 4 * N)))
    return out


def test_every_result_is_true_modulo_each_prime_or_refused():
    checked = refused = 0
    for p, r, N, a, P, k in _cases(20261023, 4000):
        try:
            R = curves.multiply(k, P, a, N)
        except curves.Degenerate:
            refused += 1
            continue
        checked += 1
        for f in (p, r):
            expected = _ref_multiply(k, (P[0] % f, P[1] % f), a % f, f)
            got = None if R is None else (R[0] % f, R[1] % f)
            assert got == expected, (N, a, P, k, f)
    assert checked > 1000 and refused > 0


def test_the_same_check_finds_the_old_gap(w164):
    """Positive control: the old (research) arithmetic violates the property on the same cases."""
    violations = 0
    for p, r, N, a, P, k in _cases(20261023, 4000):
        R = w164['mal'](k, (P[0], P[1], 1), a, N)
        for f in (p, r):
            expected = _ref_multiply(k, (P[0] % f, P[1] % f), a % f, f)
            if R is None or R[2] % f == 0:
                got = None
            else:
                zi = pow(R[2], -1, f)
                got = (R[0] * zi % f, R[1] * zi % f)
            violations += got != expected
    assert violations > 0
