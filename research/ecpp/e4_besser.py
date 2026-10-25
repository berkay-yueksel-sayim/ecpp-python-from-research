# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# e4_besser.py
# -*- coding: utf-8 -*-
# Own ECPP prover, step E4b: a faster ECPP step for large numbers (`schritt_best`, "best of K candidates").
# Reads: nothing from disk; imports e1_grundbausteine and e3_ecpp_schritt. Writes: nothing. No command-line arguments;
#   the function `schritt_best` is called by e5_teil1.py (and through it by e5_teil2.py).
# Why: `e3.schritt` takes the first usable candidate and factors m by trial division up to 2*10^5 plus Pollard rho. At 112 digits
#   q shrank by only a factor of about 2 per step (98 -> 98 digits, 39 steps in 5 min); at 477 digits the first step alone took
#   more than 6 min (rho almost always fails and costs seconds per call).
# What changed (standard knowledge, no foreign code): (1) the smooth part of m is found with gcd(m, primorial(B)) instead of a
#   loop and rho: a single gcd of a number of about 10^6 bits against m takes milliseconds; the smooth highest power comes from
#   repeated gcd(m/g, g). (2) "best of K": first collect up to K candidates (D, m, s, q), with q = m/s prime and q > bound,
#   then choose the one with the SMALLEST q (largest progress per step) and build the curve (root, point) only for that one;
#   if that fails, take the next.
# The decision stays with w164: every step is verified there. Expected outcomes:
#   F1  Chains for nextprime(10^40) and the 9 random numbers of e3 are still accepted by w164.
#   F2  Number 2 of part I (112 digits; more than 5 min and more than 39 steps with e3) is done in under 60 s.
#   F3  Progress per step improves: the median digit reduction per step at 112 digits is at least 4 (before: about 1).
import sys, math, time, random, pathlib
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent   # `HIER` = folder of this script
sys.path.insert(0, str(HIER))
import e3_ecpp_schritt as e3
import e1_grundbausteine as e1
from e1_grundbausteine import jacobi, cornacchia4, nullstelle
# gmpy2 (version 2.3) speeds up modular exponentiation about fourfold (measured 54.7 -> 12.9 ms at 1290 bits). The built-in pow
#   is replaced by a wrapper in the modules e1, e3 and this one; result and exception (ValueError for a non-invertible base)
#   stay the same. The checker w164 uses its own built-in pow and is not affected.
try:
    import gmpy2, builtins
    def _pow(a, e, n=None):
        if n is None: return builtins.pow(a, e)
        try: return int(gmpy2.powmod(a, e, n))
        except ZeroDivisionError: raise ValueError('base is not invertible for the given modulus')
    for _m in (e1, e3, sys.modules[__name__]): _m.pow = _pow
    GMPY2 = True
except ImportError:
    GMPY2 = False

B_GLATT = 3_000_000   # `B_GLATT` = smoothness bound B for the gcd-based smooth part
# `_primorial(b)`: product of all primes up to b (sieve, then a balanced product tree); `teile` = parts of the product
def _primorial(b):
    s = bytearray([1]) * (b + 1); s[0:2] = b'\x00\x00'
    for i in range(2, math.isqrt(b) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    ps = [i for i in range(b + 1) if s[i]]
    teile = [math.prod(ps[i:i + 400]) for i in range(0, len(ps), 400)]
    while len(teile) > 1:
        teile = [math.prod(teile[i:i + 2]) for i in range(0, len(teile), 2)]
    return teile[0]
PRIMORIAL = _primorial(B_GLATT)

# `glatt` = smooth: splits m into its B-smooth part (all prime powers of primes up to B) and the rest
def glatt(m):
    """(s, rest): s = B-smooth part of m with multiplicity, rest = m / s (B-rough)."""
    g = math.gcd(m, PRIMORIAL); s = 1; rest = m
    while g > 1:
        s *= g; rest //= g
        g = math.gcd(rest, g)
    return s, rest

# `schritt_best`: one ECPP step for N (probable prime). `DL` = list of discriminants D, `K` = number of candidates to collect,
# `versuche_punkt` = attempts to find a point on the curve, `verboten` = values of q that must not be chosen. Returns
# ([N, t, s, a, [x, y]], q, info) with q the next (smaller) number of the chain, or None; `info` = dict with D, A, B, m, h (class
# number) and `s_stellen` (digits of s); `rest` = cofactor q = m/s; `e3.wprim` = probable-prime test with random bases;
# `e3.pruefe_schritt` = the step checker taken from w164 (returns q if the step is valid).
def schritt_best(N, rng, DL, K=6, versuche_punkt=12, verboten=()):
    schranke = (math.isqrt(math.isqrt(N)) + 2) ** 2   # `schranke` = bound for q, about (N^(1/4) + 1)^2
    kand = []                                     # `kand` = candidates (q, D, A, B, m, s) with 4N = A^2 + D*B^2, m = group order
    for D in DL:
        if jacobi(-D % N, N) != 1: continue
        darst = cornacchia4(N, D)                 # `darst` = representation 4N = A^2 + D*B^2 (or None)
        if darst is None: continue
        A, B = darst
        for m in (N + 1 - A, N + 1 + A):
            s, rest = glatt(m)
            if rest > schranke and rest not in verboten and e3.wprim(rest, rng): kand.append((rest, D, A, B, m, s))
        if len(kand) >= K: break
    kand.sort()                                   # smallest q first
    for q, D, A, B, m, s in kand:
        f = [c % N for c in e3.hd(D)]             # Hilbert class polynomial of D modulo N (`hd`)
        j = nullstelle(f, N, rng)
        if j is None or j == 0 or j == 1728 % N: continue
        k = j * pow((1728 - j) % N, -1, N) % N    # curve y^2 = x^3 + a0*x + b0 with a0 = 3k, b0 = 2k, k = j/(1728 - j)
        a0, b0 = 3 * k % N, 2 * k % N
        for _ in range(versuche_punkt):
            x0 = rng.randrange(1, N)
            # `lam` = x0^3 + a0*x0 + b0; scaling by it puts (lam*x0, lam^2) on the curve with a = a0*lam^2, b = b0*lam^3
            lam = (x0 ** 3 + a0 * x0 + b0) % N
            if lam == 0: continue
            a, b = a0 * lam * lam % N, b0 * lam ** 3 % N
            P = (lam * x0 % N, lam * lam % N)
            t = N + 1 - m                         # `t` = trace of the curve, m = N + 1 - t
            r = e3.pruefe_schritt(N, t, s, a, P)
            if isinstance(r, int) and r == q:
                return [N, t, s, a, [P[0], P[1]]], q, dict(D=D, A=A, B=B, m=m, h=len(e3.hd(D)) - 1, s_stellen=len(str(s)))
    return None
