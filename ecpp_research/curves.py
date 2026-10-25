# SPDX-License-Identifier: Apache-2.0
"""Elliptic curve arithmetic modulo N that is correct modulo every prime divisor of N at once, or refuses.

The curve is y^2 = x^3 + a*x + b; b never enters the formulas. N need not be prime. A point is either the point at
infinity (``INFINITY``) or an affine pair (x, y) modulo N. Every case distinction of the group law is decided only when it
holds *uniformly*: an element is either 0 modulo N (hence 0 modulo every prime divisor p of N) or a unit modulo N (hence
nonzero modulo every p). If an element is 0 modulo some but not all prime divisors, the operation raises ``Degenerate``.

Invariant: every point returned is the true result modulo every prime divisor p of N (reduction modulo p is a ring
homomorphism, and every inverse taken modulo N is a unit modulo p). So a computation either gives the true multiple
modulo every p, or it raises. This is what Atkin-Morain 1993, Theorem 5.2 needs: s*P not the point at infinity, and
m*P the point at infinity, modulo every prime divisor of N.

Before this release (in the first draft of this package, and in the research checker
``research/skripte/w164_ecpp_pruefer_2026-10-02.py``) the arithmetic decided these cases modulo N only and treated degenerate points as the point at infinity; that let it
accept certificates for some composite numbers. The tests contain such certificates.
"""

from math import gcd

INFINITY = None  # the point at infinity, modulo every prime divisor of N


class Degenerate(ArithmeticError):
    """A case distinction holds modulo some prime divisors of N but not all: no group operation is defined here."""


def _unit_or_zero(d, N):
    """'zero' if d = 0 (mod N), 'unit' if gcd(d, N) = 1; raise Degenerate otherwise."""
    d %= N
    if d == 0:
        return 'zero'
    if gcd(d, N) == 1:
        return 'unit'
    raise Degenerate('a quantity vanishes modulo some prime divisors of N but not all')


def is_infinity(P, N=None):
    """True for the point at infinity."""
    return P is None


def double(P, a, N):
    """2P on y^2 = x^3 + a*x + b modulo N."""
    if P is None:
        return INFINITY
    x, y = P
    if _unit_or_zero(2 * y, N) == 'zero':  # y = 0 modulo every p: P has order 2 modulo every p
        return INFINITY
    lam = (3 * x * x + a) * pow(2 * y, -1, N) % N
    x3 = (lam * lam - 2 * x) % N
    return (x3, (lam * (x - x3) - y) % N)


def add(P, Q, a, N):
    """P + Q on y^2 = x^3 + a*x + b modulo N."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if _unit_or_zero(x2 - x1, N) == 'zero':  # equal x modulo every p
        if (y1 + y2) % N == 0:                # Q = -P modulo every p (this includes y1 = y2 = 0)
            return INFINITY
        if (y1 - y2) % N == 0:                # Q = P modulo every p
            return double(P, a, N)
        raise Degenerate('equal x, but Q = P modulo some prime divisors and Q = -P modulo others')
    lam = (y2 - y1) * pow(x2 - x1, -1, N) % N
    x3 = (lam * lam - x1 - x2) % N
    return (x3, (lam * (x1 - x3) - y1) % N)


def multiply(k, P, a, N):
    """k*P by double-and-add over the binary digits of k (k >= 0); raises Degenerate like add and double."""
    R = INFINITY
    for bit in bin(k)[2:]:
        R = double(R, a, N)
        if bit == '1':
            R = add(R, P, a, N)
    return R
