# SPDX-License-Identifier: Apache-2.0
"""Polynomials over Z/pZ and their roots (research version: ``research/ecpp/e1_grundbausteine.py``).

Polynomials are lists of coefficients, lowest first. A root found is checked (f(r) = 0 mod p) before it is returned.
"""
from .arith import inverse


def normalize(f):
    while f and f[-1] == 0:
        f.pop()
    return f


def degree(f):
    return len(f) - 1


def sub(f, g, p):
    n = max(len(f), len(g))
    return normalize([((f[i] if i < len(f) else 0) - (g[i] if i < len(g) else 0)) % p for i in range(n)])


def mul(f, g, p):
    """f*g mod p by Kronecker substitution: both polynomials as one big integer, one product, coefficients read back."""
    if not f or not g:
        return []
    bits = 2 * p.bit_length() + max(len(f), len(g)).bit_length() + 1
    F = sum(c << (i * bits) for i, c in enumerate(f))
    G = sum(c << (i * bits) for i, c in enumerate(g))
    H, mask, h = F * G, (1 << bits) - 1, []
    for _ in range(len(f) + len(g) - 1):
        h.append((H & mask) % p)
        H >>= bits
    return normalize(h)


def divmod_poly(f, g, p):
    """(quotient, remainder) of f by g mod p; ValueError if the leading coefficient of g is not invertible."""
    g = normalize(list(g))
    if not g:
        raise ZeroDivisionError
    inv = inverse(g[-1], p)
    f = list(f)
    q = [0] * max(len(f) - len(g) + 1, 1)
    for i in range(len(f) - len(g), -1, -1):
        c = f[i + len(g) - 1] * inv % p
        q[i] = c
        if c:
            for j, gj in enumerate(g):
                f[i + j] = (f[i + j] - c * gj) % p
    return normalize(q), normalize(f[:len(g) - 1])


def mod(f, g, p):
    return divmod_poly(f, g, p)[1]


def monic(f, p):
    inv = inverse(f[-1], p)
    return [c * inv % p for c in f]


def gcd_poly(f, g, p):
    f, g = normalize(list(f)), normalize(list(g))
    while g:
        f, g = g, mod(f, g, p)
    return monic(f, p) if f else f


def powmod_poly(b, e, f, p):
    """b^e modulo the polynomial f (square-and-multiply)."""
    r, b = [1], mod(b, f, p)
    for bit in bin(e)[2:]:
        r = mod(mul(r, r, p), f, p)
        if bit == '1':
            r = mod(mul(r, b, p), f, p)
    return r


def evaluate(f, x, p):
    s = 0
    for c in reversed(f):
        s = (s * x + c) % p
    return s


def linear_part(f, p):
    """Product of the distinct linear factors of f: gcd(f, x^p - x)."""
    f = monic(normalize(list(f)), p)
    xp = powmod_poly([0, 1], p, f, p)
    return gcd_poly(f, sub(xp, [0, 1], p), p)


def root(f, p, rng):
    """A root of f modulo p, or None (random equal-degree splitting with gcd(g, (x + d)^((p-1)/2) - 1))."""
    try:
        g = linear_part(f, p)
        if degree(g) < 1:
            return None
        tries = 0
        while degree(g) > 1:
            tries += 1
            if tries > 200:
                return None
            d = rng.randrange(p)
            h = powmod_poly([d, 1], (p - 1) // 2, g, p)
            k = gcd_poly(g, sub(h, [1], p), p)
            if 0 < degree(k) < degree(g):
                rest = divmod_poly(g, k, p)[0]
                g = k if degree(k) <= degree(rest) else monic(rest, p)
        r = (-g[0]) % p
    except ValueError:  # an inverse is missing: p is composite
        return None
    return r if evaluate(f, r, p) == 0 else None


__all__ = ['degree', 'evaluate', 'gcd_poly', 'linear_part', 'mul', 'powmod_poly', 'root', 'sub']
