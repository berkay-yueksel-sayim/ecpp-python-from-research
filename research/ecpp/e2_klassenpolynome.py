# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# e2_klassenpolynome.py
# -*- coding: utf-8 -*-
# Own ECPP prover, step E2: Hilbert class polynomials H_D (integer coefficients) for small |D|; `klassenpolynom(D)` is used by
# e3_ecpp_schritt.py.
# Reads: nothing from disk (imports e1_grundbausteine and mpmath). Writes: ergebnisse/e2_klassenpolynome_output.txt (next to
#   this script; the printed control report). Usage: python e2_klassenpolynome.py (no arguments; exit code 0 if all controls
#   pass).
# Basis: Atkin-Morain 1993, Sec. 2.1 (reduced forms, Prop. 2.1), Sec. 3.3 Thm 3.1 (H_D = prod (X - j(tau_C)),
#   tau_C = (-b + i*sqrt(D))/(2a)), Prop. 7.1 (H_D(0) is a cube). j is computed from standard q-expansions in TWO ways:
#   Route A: j = E4^3/Delta, E4 = 1 + 240 sum sigma_3(n) q^n, Delta = q prod (1 - q^n)^24 (product multiplied out directly);
#   Route B: Weber's f2, f2^24 = 2^12 q P(q^2)^24 / P(q)^24 with P(q) = prod (1 - q^n) as Euler's pentagonal-number series,
#            j = (f2^24 + 16)^3 / f2^24 (Atkin-Morain (7), (10), (13)).
#   The pseudocode of the paper was not copied; no foreign source code was read.
# Library: mpmath 1.3.0 (BSD license according to the package metadata); only used, its source was not read.
# Safety: H_D need NOT be proved (Atkin-Morain p. 41): a wrong H_D yields at most a curve that w164 rejects. Still, every H_D is
#   discarded if the two routes disagree, a coefficient is not close to an integer, an imaginary part does not vanish, or H_D(0)
#   is not a cube.
# Controls (the known values are the expectation):
#   K1  Class number (number of reduced forms) and genus number g = 2^(t-1) against Atkin-Morain Table 1 (p. 40), fundamental D:
#       (h, g) = (1, 1): Dmin 3, Dmax 163, count 9 | (2, 2): 15, 427, 18 | (3, 1): 23, 907, 16 | (4, 2): 39, 1555, 30 |
#       (4, 4): 84, 1435, 24 | (5, 1): 47, 2683, 25. Counted for D < 20 000; Table 1 goes up to 10^6, but its own Dmax is at most
#       2683 for these rows, so according to Table 1 none lies between 20 000 and 10^6 and the counts must be equal.
#   K2  Class number 1: H_D(X) = X - j, j = u^3 with u from Atkin-Morain p. 59 (table for Thm 8.2): D = 7: -3*5 | 11: -2^5 |
#       19: -2^5*3 | 43: -2^6*3*5 | 67: -2^5*3*5*11 | 163: -2^6*3*5*23*29; also p. 58: j(sqrt(-2)) = 20^3 (D = 8); D = 3: j = 0
#       and D = 4: j = 1728 (p. 58, curves y^2 = x^3 + b and y^2 = x^3 + ax).
#   K3  H_23(X) = X^3 + 3491750 X^2 - 5151296875 X + 23375^3 (Atkin-Morain p. 43), coefficient by coefficient.
#   K4  All fundamental D < 2000 with h <= 10: both routes agree, rounding distance < 10^-10, imaginary parts < 10^-10, H_D(0) is
#       a cube, and modulo a prime p with 4p = A^2 + D*B^2 the polynomial H_D splits into h distinct linear factors (Thm 3.2;
#       building blocks from E1).
import sys, math, hashlib, pathlib, time
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent   # `HIER` = folder of this script
sys.path.insert(0, str(HIER))
import mpmath as mp
from e1_grundbausteine import cornacchia4, linearteil, grad
ERG = HIER / 'ergebnisse'

# ---------- reduced forms ----------
# `reduzierte_formen(D, abbruch)` = reduced forms (a, b, c) of discriminant -D (`formen` = forms); `abbruch` = early exit
# as soon as more than that many forms are found.
# reduced: |b| <= a <= c, and b >= 0 if |b| = a or a = c (Prop. 2.1); a <= c and b^2 <= a^2 give D = 4ac - b^2 >= 3a^2.
def reduzierte_formen(D, abbruch=None):
    formen, a = [], 1
    while 3 * a * a <= D:
        for b in range(-a + 1, a + 1):
            if (b + D) % 2: continue                    # b = D (mod 2), otherwise b^2 + D is not divisible by 4
            z = b * b + D
            if z % (4 * a): continue
            c = z // (4 * a)
            if c < a or (a == c and b < 0): continue
            if math.gcd(math.gcd(a, abs(b)), c) != 1: continue
            formen.append((a, b, c))
            if abbruch is not None and len(formen) > abbruch: return formen
        a += 1
    return formen

def quadratfrei(n):                                   # `quadratfrei` = squarefree test by trial division
    d = 2
    while d * d <= n:
        if n % (d * d) == 0: return False
        d += 1
    return True
def fundamental(D):                                   # -D is a fundamental discriminant
    if D % 4 == 3: return quadratfrei(D)
    if D % 4 == 0: return (D // 4) % 4 in (1, 2) and quadratfrei(D // 4)
    return False
def geschlechter(D):                                  # `geschlechter` = genus number 2^(t-1), t = number of distinct primes of D
    t, n, p = 0, D, 2
    while p * p <= n:
        if n % p == 0:
            t += 1
            while n % p == 0: n //= p
        p += 1
    if n > 1: t += 1
    return 2 ** (t - 1)

# ---------- j by two routes ----------
# `n_terme` = number of series terms, `tau` = point in the upper half plane, `delta` = the discriminant function Delta
def j_route_a(tau, n_terme):
    q = mp.exp(2j * mp.pi * tau)
    e4, qn = mp.mpc(1), mp.mpc(1)
    for n in range(1, n_terme + 1):
        qn *= q
        s3 = sum(d ** 3 for d in range(1, n + 1) if n % d == 0)
        e4 += 240 * s3 * qn
    prod, qn = mp.mpc(1), mp.mpc(1)
    for n in range(1, n_terme + 1):
        qn *= q
        prod *= (1 - qn)
    delta = q * prod ** 24
    return e4 ** 3 / delta
# `fuenfeck` = pentagonal-number series (Euler)
def fuenfeck(q, n_terme):                             # P(q) = 1 + sum_{n>=1} (-1)^n (q^{n(3n-1)/2} + q^{n(3n+1)/2})
    s, n = mp.mpc(1), 1
    while n * (3 * n - 1) // 2 <= n_terme:
        v = -1 if n % 2 else 1
        s += v * (q ** (n * (3 * n - 1) // 2) + q ** (n * (3 * n + 1) // 2))
        n += 1
    return s
def j_route_b(tau, n_terme):
    q = mp.exp(2j * mp.pi * tau)
    f2_24 = 2 ** 12 * q * fuenfeck(q * q, n_terme) ** 24 / fuenfeck(q, n_terme) ** 24
    return (f2_24 + 16) ** 3 / f2_24

# ---------- H_D ----------
# `klassenpolynom(D)` returns (`koeff`, h, `qual`): integer coefficients of H_D (lowest first), class number h, and quality
# figures `qual`: `abst` = largest distance of a coefficient to the nearest integer, `imag` = largest imaginary part, `diff` =
# largest relative difference of the two routes, `dps` = working precision in digits (`stellen` = digits needed,
# `poly_a`/`poly_b` = polynomials of routes A and B).
def klassenpolynom(D):
    formen = reduzierte_formen(D)
    h = len(formen)
    stellen = sum(math.pi * math.sqrt(D) / (a * math.log(10)) for a, b, c in formen) + h * math.log10(2) + 30
    mp.mp.dps = int(stellen) + 20
    a_max = max(a for a, b, c in formen)
    n_terme = int(mp.mp.dps * math.log(10) / (2 * math.pi * math.sqrt(D) / (2 * a_max))) + 10   # |q|^N < 10^-dps
    poly_a, poly_b = [mp.mpc(1)], [mp.mpc(1)]
    for a, b, c in formen:
        tau = (-b + 1j * mp.sqrt(D)) / (2 * a)
        for poly, j in ((poly_a, j_route_a(tau, n_terme)), (poly_b, j_route_b(tau, n_terme))):
            neu = [mp.mpc(0)] * (len(poly) + 1)
            for i, k in enumerate(poly):               # polynomial times (X - j), lowest coefficient first
                neu[i + 1] += k
                neu[i] -= j * k
            poly[:] = neu
    # `neu` = new coefficient list, `skala` = scale for the relative difference of the two routes
    koeff, abst, imag, diff = [], 0.0, 0.0, 0.0
    for ka, kb in zip(poly_a, poly_b):
        r = int(mp.nint(ka.real))
        koeff.append(r)
        abst = max(abst, float(abs(ka.real - r)), float(abs(kb.real - r)))
        imag = max(imag, float(abs(ka.imag)), float(abs(kb.imag)))
        skala = max(1, abs(r))
        diff = max(diff, float(abs(ka - kb) / skala))
    return koeff, h, dict(abst=abst, imag=imag, diff=diff, dps=mp.mp.dps)

def ist_kubik(n):                                     # is |n| a perfect cube? pure integer bisection, works beyond float range
    lo, hi = 0, 1 << (abs(n).bit_length() // 3 + 2)
    while lo < hi:
        m = (lo + hi) // 2
        if m ** 3 < abs(n): lo = m + 1
        else: hi = m
    return lo ** 3 == abs(n)

def prim_klein(n):                                    # `prim_klein` = primality of a small n by trial division
    if n < 2: return False
    for d in range(2, math.isqrt(n) + 1):
        if n % d == 0: return False
    return True

# ================= controls =================
aus = []                                              # `aus` = output lines (also written to the output file)
def sag(s=''):                                        # `sag` = say: print a line and record it
    print(s, flush=True); aus.append(s)

def main():
    t0 = time.time(); gruen = {}                       # `gruen` = green: control name -> passed (True/False)
    sag('E2 — Klassenpolynome H_D · Kontrollen K1–K4 (Erwartungen im Kopf der Datei)')

    # K1: `tab1` = Table 1 rows (h, g) -> (Dmin, Dmax, count); `zaehl` = same figures counted here; `abw` = deviations
    tab1 = {(1, 1): (3, 163, 9), (2, 2): (15, 427, 18), (3, 1): (23, 907, 16), (4, 2): (39, 1555, 30), (4, 4): (84, 1435, 24),
            (5, 1): (47, 2683, 25)}
    zaehl = {}
    for D in range(3, 20000):
        if not fundamental(D): continue
        f = reduzierte_formen(D, abbruch=5)
        h = len(f)
        if h > 5: continue
        key = (h, geschlechter(D))
        mn, mx, k = zaehl.get(key, (10**9, 0, 0))
        zaehl[key] = (min(mn, D), max(mx, D), k + 1)
    abw = {k: (zaehl.get(k), v) for k, v in tab1.items() if zaehl.get(k) != v}
    gruen['K1'] = not abw
    sag(f'K1 Tab. 1 (S. 40), fundamentale D < 20 000, h ≤ 5: {len(tab1)} Zeilen verglichen, Abweichungen: {abw or "keine"} '
        f'{"✅" if gruen["K1"] else "❌"} · {time.time() - t0:.1f} s')

    # K2: `u` = cube roots of j from the table (Atkin-Morain p. 59); `soll` = expected j per D
    u = {7: -3 * 5, 11: -2**5, 19: -2**5 * 3, 43: -2**6 * 3 * 5, 67: -2**5 * 3 * 5 * 11, 163: -2**6 * 3 * 5 * 23 * 29}
    soll = {3: 0, 4: 1728, 8: 20**3}
    soll.update({D: v ** 3 for D, v in u.items()})
    fehler = []
    for D, j in sorted(soll.items()):
        koeff, h, q = klassenpolynom(D)
        if h != 1 or koeff != [-j, 1] or q['abst'] > 1e-10 or q['diff'] > 1e-10: fehler.append((D, koeff, q))
    gruen['K2'] = not fehler
    sag(f'K2 Klassenzahl 1 (S. 58–59): 9 Diskriminanten, Abweichungen: {fehler or "keine"} {"✅" if gruen["K2"] else "❌"}')

    # K3: `soll23` = expected coefficients of H_23, lowest first
    koeff, h, q = klassenpolynom(23)
    soll23 = [23375**3, -5151296875, 3491750, 1]
    gruen['K3'] = koeff == soll23 and h == 3
    sag(f'K3 H₂₃ (S. 43): berechnet {koeff[::-1]} (hoechster zuerst) · {"✅ gleich" if gruen["K3"] else "❌"} · Rundung {q["abst"]:.1e}, '
        f'Routen {q["diff"]:.1e}')

    # K4: `geprueft` = number of polynomials checked, `fehler` = failures, `gut` = numerical checks passed,
    # `zerf` = H_D splits completely modulo p
    geprueft, fehler, h_max, t1 = 0, [], 0, time.time()
    for D in range(3, 2000):
        if not fundamental(D) or len(reduzierte_formen(D, abbruch=10)) > 10: continue
        koeff, h, q = klassenpolynom(D)
        geprueft += 1; h_max = max(h_max, h)
        gut = q['abst'] < 1e-10 and q['imag'] < 1e-10 and q['diff'] < 1e-10 and ist_kubik(koeff[0])
        p = next(p for p in range(10**6 + 1, 10**7, 2) if prim_klein(p) and cornacchia4(p, D) is not None)
        zerf = grad(linearteil([k % p for k in koeff], p)) == h
        if not (gut and zerf): fehler.append((D, h, q, ist_kubik(koeff[0]), zerf))
    gruen['K4'] = not fehler and geprueft > 0
    sag(f'K4 alle fundamentalen D < 2000 mit h ≤ 10: {geprueft} Polynome (h bis {h_max}), Fehler: {fehler or "keine"} '
        f'{"✅" if gruen["K4"] else "❌"} · {time.time() - t1:.1f} s')

    alle = all(gruen.values())   # `alle` = all controls passed
    sag(f'ERGEBNIS: {"✅ alle Kontrollen gruen" if alle else "❌ " + str([k for k, v in gruen.items() if not v])} · {time.time() - t0:.1f} s')
    ERG.mkdir(exist_ok=True)
    text = '\n'.join(aus) + '\n'
    (ERG / 'e2_klassenpolynome_output.txt').write_text(text, encoding='utf-8', newline='\n')
    sag(f'md5 Ausgabe {hashlib.md5(text.encode("utf-8")).hexdigest()[:8]} · md5 Skript {hashlib.md5(pathlib.Path(__file__).read_bytes()).hexdigest()[:8]}')
    return 0 if alle else 1

if __name__ == '__main__':
    sys.exit(main())
