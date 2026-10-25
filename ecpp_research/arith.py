# SPDX-License-Identifier: Apache-2.0
"""Number-theoretic building blocks of the prover (research version: ``research/ecpp/e1_grundbausteine.py``,
``e3_ecpp_schritt.py``, ``e4_besser.py``).

None of these functions carries a proof. Every result is recomputed where it is used, and only the checker in
``verify`` decides whether a certificate is valid; an error here costs time, never correctness.
"""
import math

try:  # optional: gmpy2 speeds up modular exponentiation about fourfold at 1300 bits (research measurement)
    import gmpy2 as _gmpy2
except ImportError:  # pragma: no cover - depends on the environment
    _gmpy2 = None
HAVE_GMPY2 = _gmpy2 is not None


def powmod(a, e, n):
    """pow(a, e, n), through gmpy2 when available; raises ValueError for a non-invertible base with e < 0."""
    if _gmpy2 is None:
        return pow(a, e, n)
    try:
        return int(_gmpy2.powmod(a, e, n))
    except ZeroDivisionError:
        raise ValueError('base is not invertible for the given modulus') from None


def inverse(a, n):
    """a^(-1) mod n; ValueError if gcd(a, n) != 1."""
    return powmod(a % n, -1, n)


def jacobi(a, n):
    """Jacobi symbol (a/n) for odd n > 0, by the reciprocity law with both supplementary laws."""
    if n <= 0 or n % 2 == 0:
        raise ValueError('n must be odd and positive')
    a %= n
    t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


def sqrt_mod(a, p):
    """A square root of a modulo the (probable) prime p, or None.

    Tonelli-Shanks in the 2-Sylow subgroup of (Z/pZ)^*. The result is checked (r^2 = a mod p) before it is returned.
    """
    a %= p
    if a == 0:
        return 0
    if jacobi(a, p) != 1:
        return None
    if p % 4 == 3:
        x = powmod(a, (p + 1) // 4, p)
    else:
        q, s = p - 1, 0
        while q % 2 == 0:
            q //= 2
            s += 1
        z = 2
        while jacobi(z, p) != -1:
            z += 1
            if z > 100000:
                return None
        m, c, x, b = s, powmod(z, q, p), powmod(a, (q + 1) // 2, p), powmod(a, q, p)
        while b != 1:
            k, t = 0, b
            while t != 1:
                t = t * t % p
                k += 1
                if k >= m:
                    return None  # only possible for composite p
            e = powmod(c, 1 << (m - k - 1), p)
            x, c, b, m = x * e % p, e * e % p, b * e * e % p, k
    return x if x * x % p == a else None


def cornacchia4(N, D):
    """(A, B) with 4N = A^2 + D*B^2 for a prime N and a discriminant -D (D = 0, 3 mod 4), or None.

    Euclidean sequence of (2N, x0) with x0^2 = -D (mod 4N), stopped at the first remainder below 2*sqrt(N);
    the representation is recomputed before it is returned.
    """
    if D <= 0 or D % 4 not in (0, 3):
        raise ValueError('-D must be a discriminant (D = 0, 3 mod 4)')
    if N % 2 == 0 or math.gcd(N, D) != 1:
        return None
    r = sqrt_mod(-D % N, N)
    if r is None:
        return None
    if (r - D) % 2:
        r = N - r
    a, b = 2 * N, r
    limit = math.isqrt(4 * N)
    while b > limit:
        a, b = b, a % b
    rest = 4 * N - b * b
    if rest < 0 or rest % D:
        return None
    B = math.isqrt(rest // D)
    if B * B * D != rest:
        return None
    return b, B


def primes_up_to(n):
    """All primes <= n (sieve of Eratosthenes)."""
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b'\x00\x00'
    for i in range(2, math.isqrt(n) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(n + 1) if sieve[i]]


_SMALL = primes_up_to(300)


def probable_prime(n, rng, rounds=24):
    """Strong probable-prime test with random bases -- only for choosing candidates; it proves nothing."""
    if n < 2:
        return False
    for p in _SMALL[:50]:
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
        x = powmod(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primorial(bound):
    """Product of all primes up to ``bound`` (balanced product tree)."""
    ps = primes_up_to(bound)
    parts = [math.prod(ps[i:i + 400]) for i in range(0, len(ps), 400)]
    while len(parts) > 1:
        parts = [math.prod(parts[i:i + 2]) for i in range(0, len(parts), 2)]
    return parts[0]


def smooth_part(m, prim):
    """(s, rest): s = the part of m made of primes dividing ``prim`` (with multiplicity), rest = m / s."""
    g = math.gcd(m, prim)
    s, rest = 1, m
    while g > 1:
        s *= g
        rest //= g
        g = math.gcd(rest, g)
    return s, rest
