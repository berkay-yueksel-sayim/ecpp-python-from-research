# SPDX-License-Identifier: Apache-2.0
"""Controls of the building blocks, carried over from the research scripts e1 (K1-K5) and e2 (K1-K4)."""
import math
import os
import random

import pytest

from ecpp_research import arith, classpoly, polynomials


def _is_prime(n):
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n) + 1))


PRIMES = [p for p in range(3, 3000) if _is_prime(p)]


def test_jacobi_equals_euler():  # e1 K1
    for p in [p for p in PRIMES if p < 600]:
        for a in range(p):
            euler = 0 if a % p == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
            assert arith.jacobi(a, p) == euler


def test_sqrt_mod():  # e1 K2
    rng = random.Random(1)
    for p in PRIMES:
        for a in range(p):
            r = arith.sqrt_mod(a, p)
            residue = a % p == 0 or pow(a, (p - 1) // 2, p) == 1
            assert (r is not None and r * r % p == a) if residue else r is None
    for p in (2 ** 127 - 1, 2 ** 255 - 19, 998244353):
        for _ in range(100):
            x = rng.randrange(1, p)
            assert arith.sqrt_mod(x * x % p, p) ** 2 % p == x * x % p
    for a in range(561):  # Carmichael number: a root may be missing, but a wrong one never appears
        r = arith.sqrt_mod(a, 561)
        assert r is None or r * r % 561 == a


def test_cornacchia_against_search():  # e1 K3
    for N in PRIMES:
        for D in (3, 4, 7, 8, 11, 15, 19, 20, 23, 24, 35, 40, 43, 47, 51, 52, 67, 71, 84, 88, 163):
            if N % D == 0 or D % N == 0:
                continue
            found = None
            for B in range(math.isqrt(4 * N // D) + 1):
                A = math.isqrt(4 * N - D * B * B)
                if A * A == 4 * N - D * B * B:
                    found = (A, B)
                    break
            got = arith.cornacchia4(N, D)
            assert (got is None) == (found is None)
            if got:
                assert got[0] ** 2 + D * got[1] ** 2 == 4 * N
    assert arith.cornacchia4(17401, 8)[0] ** 2 + 8 * arith.cornacchia4(17401, 8)[1] ** 2 == 4 * 17401  # Atkin-Morain p. 59
    assert arith.cornacchia4(107, 7)[0] ** 2 + 7 * arith.cornacchia4(107, 7)[1] ** 2 == 4 * 107        # Atkin-Morain p. 60


def test_roots():  # e1 K4
    rng = random.Random(2)
    for p in [p for p in PRIMES if 50 < p < 500][::5] + [2 ** 127 - 1]:
        n1 = next(n for n in range(2, p) if arith.jacobi(n, p) == -1)
        n2 = next(n for n in range(n1 + 1, p) if arith.jacobi(n, p) == -1)
        for k in (1, 2, 3, 5):
            rs = [rng.randrange(p) for _ in range(k)]
            f = [p - n1, 0, 1]
            for r in rs:
                f = polynomials.mul(f, [(-r) % p, 1], p)
            assert polynomials.root(f, p, rng) in rs
        assert polynomials.root(polynomials.mul([p - n1, 0, 1], [p - n2, 0, 1], p), p, rng) is None


def test_h23_splitting():  # e1 K5 (Atkin-Morain p. 43, Thm 3.2; (b), (c) are our own derivation)
    H23 = [23375 ** 3, -5151296875, 3491750, 1]
    for p in [p for p in PRIMES if p > 200]:
        roots = polynomials.degree(polynomials.linear_part([c % p for c in H23], p))
        if arith.cornacchia4(p, 23) is not None:
            assert roots == 3
        elif arith.jacobi(-23 % p, p) == 1:
            assert roots == 0
        else:
            assert roots == 1


def test_class_numbers_table1():  # e2 K1 (Atkin-Morain Table 1, p. 40; counted for D < 20000)
    table = {(1, 1): (3, 163, 9), (2, 2): (15, 427, 18), (3, 1): (23, 907, 16), (4, 2): (39, 1555, 30),
             (4, 4): (84, 1435, 24), (5, 1): (47, 2683, 25)}

    def genus(D):
        t, n, p = 0, D, 2
        while p * p <= n:
            if n % p == 0:
                t += 1
                while n % p == 0:
                    n //= p
            p += 1
        return 2 ** (t - 1 + (n > 1))

    count = {}
    for D in range(3, 20000):
        if not classpoly.is_fundamental(D):
            continue
        h = len(classpoly.reduced_forms(D, stop_after=5))
        if h > 5:
            continue
        key = (h, genus(D))
        lo, hi, k = count.get(key, (10 ** 9, 0, 0))
        count[key] = (min(lo, D), max(hi, D), k + 1)
    for key, expected in table.items():
        assert count.get(key) == expected


def test_class_number_one():  # e2 K2 (Atkin-Morain pp. 58-59)
    cube_roots = {7: -3 * 5, 11: -2 ** 5, 19: -2 ** 5 * 3, 43: -2 ** 6 * 3 * 5, 67: -2 ** 5 * 3 * 5 * 11,
                  163: -2 ** 6 * 3 * 5 * 23 * 29}
    expected = {3: 0, 4: 1728, 8: 20 ** 3, **{D: u ** 3 for D, u in cube_roots.items()}}
    for D, j in expected.items():
        assert classpoly.hilbert_class_polynomial(D) == (-j, 1)


def test_h23():  # e2 K3 (Atkin-Morain p. 43)
    assert classpoly.hilbert_class_polynomial(23) == (23375 ** 3, -5151296875, 3491750, 1)


@pytest.mark.skipif(not os.environ.get('ECPP_SLOW'), reason='about a minute: set ECPP_SLOW=1')
def test_all_class_polynomials_below_2000():  # e2 K4
    p_for = {}
    for D in range(3, 2000):
        if not classpoly.is_fundamental(D) or len(classpoly.reduced_forms(D, stop_after=10)) > 10:
            continue
        H = classpoly.hilbert_class_polynomial(D)  # raises if the two routes disagree
        p = p_for.setdefault(D, next(p for p in range(10 ** 6 + 1, 10 ** 7, 2) if _is_prime(p) and arith.cornacchia4(p, D)))
        assert polynomials.degree(polynomials.linear_part([c % p for c in H], p)) == len(H) - 1
