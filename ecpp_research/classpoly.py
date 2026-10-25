# SPDX-License-Identifier: Apache-2.0
"""Hilbert class polynomials H_D for small |D| (research version: ``research/ecpp/e2_klassenpolynome.py``).

Basis (Atkin-Morain 1993): reduced forms (Sec. 2.1, Prop. 2.1); H_D = prod (X - j(tau_C)) over the reduced forms,
tau_C = (-b + i*sqrt(D)) / (2a) (Sec. 3.3, Thm 3.1); H_D(0) is a cube (Prop. 7.1). j is computed in two ways:
route A, j = E4^3 / Delta from the q-expansions; route B, Weber's f2 with j = (f2^24 + 16)^3 / f2^24.
A polynomial is rejected unless both routes agree, every coefficient is close to an integer, the imaginary parts
vanish, and H_D(0) is a cube. H_D needs no proof (Atkin-Morain p. 41): a wrong H_D can only produce a curve that
the checker rejects.
"""
import math
from functools import lru_cache

import mpmath as mp


def reduced_forms(D, stop_after=None):
    """Reduced primitive forms (a, b, c) of discriminant -D: |b| <= a <= c, b >= 0 if |b| = a or a = c."""
    forms, a = [], 1
    while 3 * a * a <= D:
        for b in range(-a + 1, a + 1):
            if (b + D) % 2:
                continue
            z = b * b + D
            if z % (4 * a):
                continue
            c = z // (4 * a)
            if c < a or (a == c and b < 0):
                continue
            if math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            forms.append((a, b, c))
            if stop_after is not None and len(forms) > stop_after:
                return forms
        a += 1
    return forms


def class_number(D):
    return len(reduced_forms(D))


def _squarefree(n):
    d = 2
    while d * d <= n:
        if n % (d * d) == 0:
            return False
        d += 1
    return True


def is_fundamental(D):
    """-D is a fundamental discriminant."""
    if D % 4 == 3:
        return _squarefree(D)
    if D % 4 == 0:
        return (D // 4) % 4 in (1, 2) and _squarefree(D // 4)
    return False


def _j_route_a(tau, terms):
    q = mp.exp(2j * mp.pi * tau)
    e4, qn = mp.mpc(1), mp.mpc(1)
    for n in range(1, terms + 1):
        qn *= q
        e4 += 240 * sum(d ** 3 for d in range(1, n + 1) if n % d == 0) * qn
    prod, qn = mp.mpc(1), mp.mpc(1)
    for n in range(1, terms + 1):
        qn *= q
        prod *= (1 - qn)
    return e4 ** 3 / (q * prod ** 24)


def _pentagonal(q, terms):
    s, n = mp.mpc(1), 1
    while n * (3 * n - 1) // 2 <= terms:
        sign = -1 if n % 2 else 1
        s += sign * (q ** (n * (3 * n - 1) // 2) + q ** (n * (3 * n + 1) // 2))
        n += 1
    return s


def _j_route_b(tau, terms):
    q = mp.exp(2j * mp.pi * tau)
    f2_24 = 2 ** 12 * q * _pentagonal(q * q, terms) ** 24 / _pentagonal(q, terms) ** 24
    return (f2_24 + 16) ** 3 / f2_24


def _is_cube(n):
    lo, hi = 0, 1 << (abs(n).bit_length() // 3 + 2)
    while lo < hi:
        m = (lo + hi) // 2
        if m ** 3 < abs(n):
            lo = m + 1
        else:
            hi = m
    return lo ** 3 == abs(n)


class ClassPolynomialError(ArithmeticError):
    """The two routes disagree, or a coefficient is not an integer, or H_D(0) is not a cube."""


@lru_cache(maxsize=None)
def hilbert_class_polynomial(D):
    """Integer coefficients of H_D (lowest first), checked as described in the module docstring."""
    forms = reduced_forms(D)
    h = len(forms)
    digits = sum(math.pi * math.sqrt(D) / (a * math.log(10)) for a, b, c in forms) + h * math.log10(2) + 30
    with mp.workdps(int(digits) + 20):
        a_max = max(a for a, b, c in forms)
        terms = int(mp.mp.dps * math.log(10) / (2 * math.pi * math.sqrt(D) / (2 * a_max))) + 10
        poly_a, poly_b = [mp.mpc(1)], [mp.mpc(1)]
        for a, b, c in forms:
            tau = (-b + 1j * mp.sqrt(D)) / (2 * a)
            for poly, j in ((poly_a, _j_route_a(tau, terms)), (poly_b, _j_route_b(tau, terms))):
                new = [mp.mpc(0)] * (len(poly) + 1)
                for i, k in enumerate(poly):
                    new[i + 1] += k
                    new[i] -= j * k
                poly[:] = new
        coeffs, dist, imag, diff = [], 0.0, 0.0, 0.0
        for ka, kb in zip(poly_a, poly_b):
            r = int(mp.nint(ka.real))
            coeffs.append(r)
            dist = max(dist, float(abs(ka.real - r)), float(abs(kb.real - r)))
            imag = max(imag, float(abs(ka.imag)), float(abs(kb.imag)))
            diff = max(diff, float(abs(ka - kb) / max(1, abs(r))))
    if dist > 1e-10 or imag > 1e-10 or diff > 1e-10 or not _is_cube(coeffs[0]):
        raise ClassPolynomialError(f'H_{D}: rounding {dist:.1e}, imaginary {imag:.1e}, routes {diff:.1e}, '
                                   f'cube {_is_cube(coeffs[0])}')
    return tuple(coeffs)


def discriminants(h_max=8, d_max=6000):
    """Fundamental D >= 7 below d_max with class number <= h_max, sorted by (h, D); D = 3, 4 (j = 0, 1728) are left out."""
    found = []
    for D in range(7, d_max):
        if D % 4 not in (0, 3) or not is_fundamental(D):
            continue
        h = len(reduced_forms(D, stop_after=h_max))
        if h <= h_max:
            found.append((h, D))
    found.sort()
    return [D for h, D in found]
