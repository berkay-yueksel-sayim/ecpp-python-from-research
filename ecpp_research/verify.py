# SPDX-License-Identifier: Apache-2.0
"""Independent checker for ECPP certificates (pure Python, no outside library).

Theorem used (Goldwasser-Kilian 1986; Atkin-Morain 1993, Theorem 5.2 and Corollary 5.1):
let N > 1 with gcd(N, 6) = 1, let E: y^2 = x^3 + a*x + b with gcd(4a^3 + 27b^2, N) = 1, let
m = s*q with q prime and q > (N^(1/4) + 1)^2, and let P be a point such that, modulo every prime
divisor p of N, s*P is not the point at infinity and m*P is. Then N is prime.
(For p | N the point Q = s*P has order q modulo p, so q <= #E(F_p) <= (sqrt(p) + 1)^2, and
q > (N^(1/4) + 1)^2 forces p > sqrt(N).)

What is computed (see ``curves``: every operation is the true group operation modulo every p | N, or it raises
``curves.Degenerate``, which rejects the step):
  (i)  Q = s*P is an affine point, so it is not the point at infinity modulo any p | N;
  (ii) q*Q is the point at infinity, modulo every p | N.
The chain is linked by N_{i+1} = q_i; the last q is at most 2^64 and is proved prime by the
deterministic Miller-Rabin test with the first twelve prime bases. The bound check is written as
q > (isqrt(isqrt(N)) + 2)^2, which is stricter than q > (N^(1/4) + 1)^2.

The conditions follow the research checker that verified the elliptic curve certificates of the paper series
(``research/skripte/w164_ecpp_pruefer_2026-10-02.py``), in the same order. Its curve arithmetic decided case
distinctions modulo N only; that let it accept certificates for some composite numbers. This checker instead refuses
every operation that is not uniform modulo all prime divisors of N.
"""
from dataclasses import dataclass
from math import gcd, isqrt

from . import curves

_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime_small(n):
    """Deterministic Miller-Rabin test with the first twelve prime bases.

    Correct for all n < psi_12 = 318665857834031151167461, the smallest strong pseudoprime to these
    bases (Sorenson and Webster, Math. Comp. 86 (2017), Theorem 1.1); the checker uses it only for
    n <= 2^64.
    """
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def _malformed(step):
    """Return a reason if ``step`` is not (N, t, s, a, (x, y)) with ints, else None."""
    if not isinstance(step, (tuple, list)) or len(step) != 5:
        return 'malformed step (expected N, t, s, a, (x, y))'
    if not all(_is_int(v) for v in step[:4]):
        return 'malformed step (N, t, s and a must be integers)'
    P = step[4]
    if not isinstance(P, (tuple, list)) or len(P) != 2 or not all(_is_int(v) for v in P):
        return 'malformed step (the point must be a pair of integers)'
    return None


def check_step(N, t, s, a, P):
    """Check one step. Returns q = m/s (an int) if the step is valid, otherwise a str naming the failed condition."""
    bad = _malformed((N, t, s, a, P))
    if bad:
        return bad
    x, y = P
    if N < 2 or gcd(N, 6) != 1:
        return 'gcd(N, 6) != 1'
    if t * t >= 4 * N:
        return 't^2 >= 4N'
    m = N + 1 - t
    if s <= 0 or m % s:
        return 's does not divide m'
    q = m // s
    r = isqrt(isqrt(N)) + 1
    if q <= (r + 1) ** 2:
        return 'q too small'
    b = (y * y - x * x * x - a * x) % N
    if gcd((4 * a ** 3 + 27 * b * b) % N, N) != 1:
        return 'curve singular modulo a divisor of N'
    try:
        Q = curves.multiply(s, (x % N, y % N), a, N)
        if Q is curves.INFINITY:
            return 's*P is the point at infinity'
        R = curves.multiply(q, Q, a, N)
    except curves.Degenerate:
        return 'degenerate group operation modulo a divisor of N'
    if R is not curves.INFINITY:
        return 'q*Q is not the point at infinity'
    return q


@dataclass
class Result:
    """Outcome of ``check_certificate``: ``valid``, the number N proved prime, the number of steps, and the reason if invalid."""
    valid: bool
    number: int
    steps: int
    reason: str = ''

    def __bool__(self):
        return self.valid


def check_certificate(steps, number=None):
    """Check a whole certificate (list of steps (N, t, s, a, (x, y))).

    If ``number`` is given, the first N must equal it. Returns a ``Result``; input that is not a list of such steps
    gives an invalid ``Result`` with the reason, not an exception.
    """
    if not isinstance(steps, (tuple, list)):
        return Result(False, 0, 0, 'not a list of steps')
    steps = list(steps)
    if not steps:
        return Result(False, 0, 0, 'empty certificate')
    for i, step in enumerate(steps):
        bad = _malformed(step)
        if bad:
            first = steps[0][0] if i else 0
            return Result(False, first, len(steps), f'step {i}: {bad}')
    first = steps[0][0]
    if number is not None and first != number:
        return Result(False, first, len(steps), 'the certificate belongs to a different number')
    expected = None
    for i, (N, t, s, a, P) in enumerate(steps):
        if expected is not None and N != expected:
            return Result(False, first, len(steps), f'step {i}: chain broken (N != q of the previous step)')
        q = check_step(N, t, s, a, P)
        if isinstance(q, str):
            return Result(False, first, len(steps), f'step {i}: {q}')
        expected = q
    if expected > 2 ** 64 or not is_prime_small(expected):
        return Result(False, first, len(steps), 'last q is not a prime <= 2^64')
    return Result(True, first, len(steps))
