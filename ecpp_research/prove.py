# SPDX-License-Identifier: Apache-2.0
"""The prover: builds an ECPP certificate chain from N down to a prime q <= 2^64.

Research version: ``research/ecpp/e3_ecpp_schritt.py`` (one step), ``e4_besser.py`` (best of K candidates, smooth
part by a gcd with a primorial) and ``e5_teil1.py`` (chain with one-stage backtracking, discriminant lists that
grow when no step is found). Basis (Atkin-Morain 1993, written from the published description, no foreign code):
Sec. 5.3 (choose D with 4N = A^2 + D*B^2, candidates m = N + 1 -/+ A, m = s*q with q > (N^(1/4) + 1)^2),
Sec. 8.6.2 Eq. (39) (the curve y^2 = x^3 + 3kc^2 x + 2kc^3 with k = j/(1728 - j) has invariant j), Sec. 8.6.3
(a point without square roots: with lam = x0^3 + a0*x0 + b0 the point (lam*x0, lam^2) lies on
Y^2 = X^3 + a0*lam^2*X + b0*lam^3; the twist is covered by the choice of x0).

The prover only SEARCHES. Every step it returns has been accepted by ``verify.check_step``, and
``prove`` checks the whole chain with ``verify.check_certificate`` before returning it, so an error here can
never produce a false proof.
"""
import math
import random
import time

from . import arith, classpoly, polynomials
from .verify import check_certificate, check_step

SMOOTH_BOUND = 3_000_000
_PRIMORIAL = None


def _primorial():
    global _PRIMORIAL
    if _PRIMORIAL is None:
        _PRIMORIAL = arith.primorial(SMOOTH_BOUND)
    return _PRIMORIAL


def _class_polynomial_mod(D, N):
    try:
        return [c % N for c in classpoly.hilbert_class_polynomial(D)]
    except classpoly.ClassPolynomialError:
        return None


def step(N, rng, discriminants, K=6, point_tries=12, forbidden=()):
    """One ECPP step for N (a probable prime > 2^64).

    Collects up to K candidates (q, D, A, B, m, s) with m = N + 1 -/+ A, s the smooth part of m and q = m/s a probable
    prime above the bound, tries them with the smallest q first, and returns ([N, t, s, a, [x, y]], q, info) for the
    first one whose step ``verify.check_step`` accepts; None if there is none. ``forbidden`` lists values of q to skip.
    """
    bound = (math.isqrt(math.isqrt(N)) + 2) ** 2
    cand = []
    for D in discriminants:
        if arith.jacobi(-D % N, N) != 1:
            continue
        rep = arith.cornacchia4(N, D)
        if rep is None:
            continue
        A, B = rep
        for m in (N + 1 - A, N + 1 + A):
            s, rest = arith.smooth_part(m, _primorial())
            if rest > bound and rest not in forbidden and arith.probable_prime(rest, rng):
                cand.append((rest, D, A, B, m, s))
        if len(cand) >= K:
            break
    cand.sort()
    for q, D, A, B, m, s in cand:
        f = _class_polynomial_mod(D, N)
        if f is None:
            continue
        j = polynomials.root(f, N, rng)
        if j is None or j == 0 or j == 1728 % N:
            continue
        try:
            k = j * arith.inverse((1728 - j) % N, N) % N
        except ValueError:  # 1728 - j shares a factor with N: N is composite
            return None
        a0, b0 = 3 * k % N, 2 * k % N
        t = N + 1 - m
        for _ in range(point_tries):
            x0 = rng.randrange(1, N)
            lam = (x0 ** 3 + a0 * x0 + b0) % N
            if lam == 0:
                continue
            a = a0 * lam * lam % N
            P = (lam * x0 % N, lam * lam % N)
            if check_step(N, t, s, a, P) == q:
                return [N, t, s, a, [P[0], P[1]]], q, dict(D=D, h=len(f) - 1, m=m, s_digits=len(str(s)))
    return None


class ProofNotFound(RuntimeError):
    """No certificate within the limits (time, backtracking). This says nothing about N; it is a limit of the search."""


# discriminant lists of growing size, as in the research runs: (h_max, d_max)
DISCRIMINANT_LEVELS = ((8, 6000), (14, 40000), (24, 200000))


def prove(N, time_limit=3600.0, seed=None, max_backtracks=40, log=None):
    """Certificate (list of steps) proving that N > 2^64 is prime.

    Raises ValueError if N <= 2^64 (use ``is_prime_small``) or N is not a probable prime, and ProofNotFound if the
    search hits the time limit or the backtracking limit. ``log`` is an optional callable receiving progress lines.
    """
    if N <= 2 ** 64:
        raise ValueError('N <= 2^64: use is_prime_small, which is deterministic there')
    rng = random.Random(seed)
    if not arith.probable_prime(N, rng):
        raise ValueError('N is not a probable prime')
    t0, stages, forbidden, backtracks, n = time.time(), [], {}, 0, N
    lists = [None] * len(DISCRIMINANT_LEVELS)
    while n > 2 ** 64:
        if time.time() - t0 >= time_limit:
            raise ProofNotFound(f'time limit {time_limit:.0f} s reached at {len(str(n))} digits')
        r = None
        for level, (h_max, d_max) in enumerate(DISCRIMINANT_LEVELS):
            if lists[level] is None:
                lists[level] = classpoly.discriminants(h_max, d_max)
            r = step(n, rng, lists[level], forbidden=forbidden.get(n, set()))
            if r is not None:
                break
        if r is None:
            if not stages or backtracks >= max_backtracks:
                raise ProofNotFound(f'no step found at {len(str(n))} digits ({backtracks} backtracks)')
            st, q, info = stages.pop()
            backtracks += 1
            forbidden.setdefault(st[0], set()).add(q)
            n = st[0]
            if log:
                log(f'no step at {len(str(q))} digits, back to {len(str(n))} digits (backtrack {backtracks})')
            continue
        st, q, info = r
        stages.append((st, q, info))
        if log:
            log(f'step {len(stages)}: {len(str(n))} -> {len(str(q))} digits, D = {info["D"]} (h = {info["h"]}), '
                f'{time.time() - t0:.0f} s')
        n = q
    steps = [tuple(s[0][:4]) + (tuple(s[0][4]),) for s in stages]
    result = check_certificate(steps, N)
    if not result.valid:  # cannot happen if every step was checked; kept as the final guard
        raise ProofNotFound(f'assembled chain rejected: {result.reason}')
    return steps
