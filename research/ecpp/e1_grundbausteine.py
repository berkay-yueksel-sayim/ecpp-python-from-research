# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# e1_grundbausteine.py
# -*- coding: utf-8 -*-
# Own ECPP prover, E1: basic building blocks (Jacobi symbol, square root modulo p, representation 4N = A² + D·B²,
# polynomial arithmetic over Z/pZ, root finding modulo p) together with their controls.
#
# Reads:    nothing. Writes: ergebnisse/e1_grundbausteine_output.txt (relative to this file), and printed output.
# Usage:    python e1_grundbausteine.py   (no arguments; pure Python, no outside library); exit code 0 if all controls are green.
# Controls: K1 to K5 (listed below); PK = positive control, NK = negative control.
#
# BASIS: standard number theory (Euler criterion, quadratic reciprocity law, cyclic group (Z/pZ)^*, x^p − x = Π(x − a) over F_p).
#   Atkin–Morain 1993 name the building blocks in §8.4.2 (representation 4N = A² + D·B²), §8.4.3 (square root modulo N) and §8.6.1
#   (root modulo N). Their pseudocode was NOT copied; no foreign source code was read. Written from the published description.
# SAFETY: no building block here carries a proof. Every result is recomputed on the spot (r² ≡ a; 4N = A² + D·B²; f(r) ≡
#   0) and discarded otherwise. The primality proof is decided in the end solely by the checker w164. An error here
#   costs time, never correctness.
# EXPECTATIONS (set before the run; the controls are below):
#   K1  jacobi(a, p) = Euler criterion a^((p−1)/2) for all primes 3 ≤ p < 600 and all a.
#   K2  `wurzel`: for primes p < 3000 (all residues p mod 8) and every a: quadratic residue ⇒ r² ≡ a; non-residue ⇒ None. Large: p
#       = 2^127 − 1 (≡ 3 mod 4), p = 2^255 − 19 (≡ 5 mod 8), p = 998244353 = 119·2^23 + 1 (2-part 2^23, the hardest case for the
#       loop). NK: n = 561 (Carmichael number): every answer ≠ None satisfies r² ≡ a (mod n); never a wrong root.
#   K3  cornacchia4(N, D) against exhaustive search: for all primes N < 3000 with N ∤ D and all D from the list below, a solution
#       exactly if the search finds one; every solution satisfies 4N = A² + D·B². Examples from Atkin–Morain pp. 59–60:
#       17401 with D = 8 (there 17401 = 101² + 2·60², so 4N = 202² + 8·60²) and 107 with D = 7 (there 107 = 10² + 7·1²).
#   K4  `nullstelle`: f = Π(x − r_i) · (x² − n), n a non-residue ⇒ the root found is among the r_i; NK f = (x² − n)(x² − n') with
#       two non-residues ⇒ None. Small (p < 500) and large (p = 2^127 − 1).
#   K5  Cross-check with H₂₃(X) = X³ + 3491750X² − 5151296875X + 23375³ (Atkin–Morain p. 43) for primes 200 < p < 3000:
#       (a) 4p = A² + 23B² solvable ⇒ H₂₃ splits modulo p into three distinct linear factors (Theorem 3.2 there);
#       (b) (−23/p) = 1, but 4p = A² + 23B² unsolvable ⇒ no root [plausible: the Frobenius element is then a 3-cycle];
#       (c) (−23/p) = −1 ⇒ exactly one root [plausible: the Frobenius element is then a reflection, it fixes
#           exactly one of the three].
#       (b) and (c) are our own derivation, not from the paper; if they fail, the derivation is suspect first, not the code.
import sys, math, random, hashlib, pathlib, time
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent   # `HIER` = this directory
ERG = HIER / 'ergebnisse'                        # `ERG` = results directory

# ---------- Jacobi symbol (reciprocity law with both supplementary laws) ----------
def jacobi(a, n):   # `jacobi` = Jacobi symbol (a/n) for odd positive n; returns 1, -1 or 0
    if n <= 0 or n % 2 == 0: raise ValueError('n muss ungerade und positiv sein')
    a %= n; t = 1
    while a:
        while a % 2 == 0:                      # (2/n) = −1 exactly for n ≡ 3, 5 (mod 8)
            a //= 2
            if n % 8 in (3, 5): t = -t
        a, n = n, a                            # reciprocity: sign change exactly when both are ≡ 3 (mod 4)
        if a % 4 == 3 and n % 4 == 3: t = -t
        a %= n
    return t if n == 1 else 0

# ---------- Square root modulo p ----------
# Derivation: (Z/pZ)^* is cyclic of order p − 1 = 2^s·q, q odd. For a quadratic residue a, x = a^((q+1)/2) is a root of a·b
#   with b = a^q; b lies in the 2-Sylow subgroup and has order 2^k there, k < s. c = z^q (z a non-residue) generates the
#   2-Sylow subgroup. Multiplying x by e = c^(2^(m−k−1)) multiplies b by e², and the order of b strictly decreases. After at
#   most s steps b = 1, hence x² = a.
def wurzel(a, p):   # `wurzel` = square root of a modulo the (probable) prime p: returns a root, or None if there is none
    a %= p
    if a == 0: return 0
    if jacobi(a, p) != 1: return None
    if p % 4 == 3:
        x = pow(a, (p + 1) // 4, p)
    else:
        q, s = p - 1, 0
        while q % 2 == 0: q //= 2; s += 1
        z = 2
        while jacobi(z, p) != -1:
            z += 1
            if z > 100000: return None
        m, c, x, b = s, pow(z, q, p), pow(a, (q + 1) // 2, p), pow(a, q, p)
        while b != 1:
            k, t = 0, b
            while t != 1:
                t = t * t % p; k += 1
                if k >= m: return None         # only possible for composite p
            e = pow(c, 1 << (m - k - 1), p)
            x, c, b, m = x * e % p, e * e % p, b * e * e % p, k
    return x if x * x % p == a else None

# ---------- Representation 4N = A² + D·B² ----------
# Derivation (Euclid/lattice argument): if 4N = A² + D·B² with N prime, N ∤ D, then A ≡ ±x0·B (mod 2N) for a root x0 of −D
#   modulo 4N; the remainders of the Euclidean sequence of (2N, x0) run through the short vectors of this lattice, and the
#   first remainder below 2√N is the sought A. The validity is NOT used here as a theorem but checked in K3 against
#   exhaustive search; every output is recomputed.
def cornacchia4(N, D):   # `cornacchia4` = finds (A, B) with 4N = A² + D·B² for a prime N and a discriminant −D, or returns None
    if D <= 0 or D % 4 not in (0, 3): raise ValueError('−D muss eine Diskriminante sein (D ≡ 0, 3 mod 4)')
    if N % 2 == 0 or math.gcd(N, D) != 1: return None
    r = wurzel(-D % N, N)
    if r is None: return None
    if (r - D) % 2: r = N - r                  # x0 ≡ D (mod 2) ⇒ x0² ≡ −D (mod 4), and with r² ≡ −D (mod N): x0² ≡ −D (mod 4N)
    a, b = 2 * N, r
    grenze = math.isqrt(4 * N)                 # `grenze` = bound at which the Euclidean sequence stops
    while b > grenze:
        a, b = b, a % b
    rest = 4 * N - b * b
    if rest < 0 or rest % D: return None
    B = math.isqrt(rest // D)
    if B * B * D != rest: return None
    return b, B

# ---------- Polynomials over Z/pZ (lists, lowest coefficient first) ----------
def p_norm(f):   # `p_norm` = remove leading zero coefficients (in place) and return the list
    while f and f[-1] == 0: f.pop()
    return f
def grad(f): return len(f) - 1   # `grad` = degree (-1 for the zero polynomial)
def p_sub(f, g, p):   # `p_sub` = f − g modulo p
    n = max(len(f), len(g))
    return p_norm([((f[i] if i < len(f) else 0) - (g[i] if i < len(g) else 0)) % p for i in range(n)])
def p_mul(f, g, p):   # `p_mul` = product f·g modulo p
    # Kronecker substitution: both polynomials as one big integer, one product, then read the coefficients back
    if not f or not g: return []
    bits = 2 * p.bit_length() + max(len(f), len(g)).bit_length() + 1
    F = sum(c << (i * bits) for i, c in enumerate(f)); G = sum(c << (i * bits) for i, c in enumerate(g))
    H = F * G; maske = (1 << bits) - 1; h = []
    for _ in range(len(f) + len(g) - 1):
        h.append((H & maske) % p); H >>= bits
    return p_norm(h)
def p_divmod(f, g, p):   # `p_divmod` = polynomial division with remainder modulo p: returns (quotient, remainder)
    g = p_norm(list(g))
    if not g: raise ZeroDivisionError
    inv = pow(g[-1], -1, p)                    # ValueError if the leading coefficient is not invertible (composite p)
    f = list(f); q = [0] * max(len(f) - len(g) + 1, 1)
    for i in range(len(f) - len(g), -1, -1):
        c = f[i + len(g) - 1] * inv % p
        q[i] = c
        if c:
            for j, gj in enumerate(g):
                f[i + j] = (f[i + j] - c * gj) % p
    return p_norm(q), p_norm(f[:len(g) - 1])
def p_mod(f, g, p): return p_divmod(f, g, p)[1]   # `p_mod` = remainder of f modulo g
def p_monisch(f, p):   # `p_monisch` = f scaled to be monic
    inv = pow(f[-1], -1, p)
    return [c * inv % p for c in f]
def p_ggt(f, g, p):   # `p_ggt` = monic greatest common divisor (Euclidean algorithm)
    f, g = p_norm(list(f)), p_norm(list(g))
    while g:
        f, g = g, p_mod(f, g, p)
    return p_monisch(f, p) if f else f
def p_potmod(b, e, f, p):   # `p_potmod` = b^e modulo the polynomial f (square-and-multiply over the binary digits of e)
    r, b = [1], p_mod(b, f, p)
    for bit in bin(e)[2:]:
        r = p_mod(p_mul(r, r, p), f, p)
        if bit == '1': r = p_mod(p_mul(r, b, p), f, p)
    return r
def p_wert(f, x, p):   # `p_wert` = value f(x) modulo p (Horner scheme)
    s = 0
    for c in reversed(f): s = (s * x + c) % p
    return s

# ---------- Root modulo p ----------
# Derivation: x^p − x = Π_{a ∈ F_p} (x − a), hence g = gcd(f, x^p − x) is the product of the distinct linear factors of f. For a
#   root a of g and a δ with a + δ ≠ 0, (a + δ)^((p−1)/2) = ±1 depending on whether a + δ is a quadratic residue (Euler
#   criterion). gcd(g, (x + δ)^((p−1)/2) − 1) collects exactly the roots with quadratic-residue a + δ; for random δ this
#   separates two distinct roots with probability about 1/2. One keeps the smaller proper divisor until degree 1 is reached.
def linearteil(f, p):   # `linearteil` = linear part of f: the product of its distinct linear factors, gcd(f, x^p − x)
    f = p_monisch(p_norm(list(f)), p)
    xp = p_potmod([0, 1], p, f, p)
    return p_ggt(f, p_sub(xp, [0, 1], p), p)
# `nullstelle` = a root of f modulo p (random splitting with the generator `rng`), or None if there is none
def nullstelle(f, p, rng):
    try:
        g = linearteil(f, p)
        if grad(g) < 1: return None
        versuche = 0   # `versuche` = number of splitting attempts so far
        while grad(g) > 1:
            versuche += 1
            if versuche > 200: return None
            d = rng.randrange(p)
            h = p_potmod([d, 1], (p - 1) // 2, g, p)
            k = p_ggt(g, p_sub(h, [1], p), p)
            if 0 < grad(k) < grad(g):
                rest = p_divmod(g, k, p)[0]
                g = k if grad(k) <= grad(rest) else p_monisch(rest, p)
        r = (-g[0]) % p
    except ValueError:                         # an inverse is missing: p is composite
        return None
    return r if p_wert(f, r, p) == 0 else None

# ================= Controls =================
aus = []   # `aus` = lines of the output file
def sag(s=''):   # `sag` = print a line and record it for the output file
    print(s, flush=True); aus.append(s)
def ist_prim(n):   # primality test by trial division (only used for small n)
    if n < 2: return False
    for d in range(2, math.isqrt(n) + 1):
        if n % d == 0: return False
    return True
def main():
    # `gruen` = status ("green" = passed) of each control; `PR` = the primes 3 <= p < 3000
    t0 = time.time(); rng = random.Random(20261003); gruen = {}
    PR = [p for p in range(3, 3000) if ist_prim(p)]
    sag('E1 — Grundbausteine des eigenen ECPP-Beweisers · Kontrollen K1–K5 (Erwartungen im Kopf der Datei)')

    ok = all(jacobi(a, p) == (0 if a % p == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)) for p in PR if p < 600 for a in range(p))
    gruen['K1'] = ok; sag(f'K1 Jacobi = Euler, alle p < 600, alle a: {"✅" if ok else "❌"}')

    # K2: `fehler` = errors for the small primes, `anzahl` = number of cases, `gfehler` = errors for the large
    #     primes, `nk` = NK for n = 561
    fehler, anzahl = 0, 0
    for p in PR:
        for a in range(p):
            r = wurzel(a, p); anzahl += 1
            qr = a % p == 0 or pow(a, (p - 1) // 2, p) == 1
            if qr and (r is None or r * r % p != a % p): fehler += 1
            if not qr and r is not None: fehler += 1
    gross = [(2**127 - 1, '2^127−1'), (2**255 - 19, '2^255−19'), (998244353, '119·2^23+1')]
    gfehler = 0
    for p, name in gross:
        for _ in range(200):
            x = rng.randrange(1, p); a = x * x % p
            r = wurzel(a, p)
            if r is None or r * r % p != a: gfehler += 1
            n = rng.randrange(1, p)
            if jacobi(n, p) == -1 and wurzel(n, p) is not None: gfehler += 1
    nk = all((lambda r: r is None or r * r % 561 == a)(wurzel(a, 561)) for a in range(561))
    ok = fehler == 0 and gfehler == 0 and nk
    gruen['K2'] = ok
    sag(f'K2 Wurzel: {anzahl} Faelle p < 3000 mit {fehler} Fehlern · gross (2^127−1, 2^255−19, 119·2^23+1, je 200 + 200 NK): {gfehler} Fehler · '
        f'NK 561 nie falsche Wurzel: {"✅" if nk else "❌"} ⇒ {"✅" if ok else "❌"}')

    # K3: `DLISTE` = list of discriminants D; `mit` / `ohne` = cases with / without a solution; `fehler` =
    #     deviations from the search
    DLISTE = [3, 4, 7, 8, 11, 15, 19, 20, 23, 24, 35, 40, 43, 47, 51, 52, 67, 71, 84, 88, 163]
    fehler, mit, ohne = 0, 0, 0
    for N in PR:
        for D in DLISTE:
            if N % D == 0 or D % N == 0: continue
            such = None   # `such` = result (A, B) of the exhaustive search, or None
            for B in range(0, math.isqrt(4 * N // D) + 1):
                A2 = 4 * N - D * B * B; A = math.isqrt(A2)
                if A * A == A2: such = (A, B); break
            erg = cornacchia4(N, D)
            if (erg is None) != (such is None): fehler += 1
            if erg is not None:
                mit += 1
                if erg[0] ** 2 + D * erg[1] ** 2 != 4 * N: fehler += 1
            else: ohne += 1
    b1 = cornacchia4(17401, 8); b2 = cornacchia4(107, 7)   # the two examples from Atkin–Morain pp. 59–60
    bsp = b1 is not None and b1[0] ** 2 + 8 * b1[1] ** 2 == 4 * 17401 and b2 is not None and b2[0] ** 2 + 7 * b2[1] ** 2 == 4 * 107
    ok = fehler == 0 and bsp and mit > 0 and ohne > 0
    gruen['K3'] = ok
    sag(f'K3 Darstellung 4N = A² + D·B² gegen erschoepfende Suche: {mit} mit, {ohne} ohne Loesung, {fehler} Abweichungen · '
        f'Atkin–Morain S. 59–60: 17401/D=8 → {b1}, 107/D=7 → {b2} ⇒ {"✅" if ok else "❌"}')

    # K4: `faelle` = number of cases; `n1`, `n2` = two non-residues; `rs` = random roots r_i; `f` = Π(x − r_i)·(x² − n1)
    fehler, faelle = 0, 0
    for p in [p for p in PR if 50 < p < 500][::5] + [2**127 - 1]:
        n1 = next(n for n in range(2, p) if jacobi(n, p) == -1)
        n2 = next(n for n in range(n1 + 1, p) if jacobi(n, p) == -1)
        for k in (1, 2, 3, 5):
            rs = [rng.randrange(p) for _ in range(k)]
            f = [p - n1, 0, 1]
            for r in rs: f = p_mul(f, [(-r) % p, 1], p)
            z = nullstelle(f, p, rng); faelle += 1
            if z is None or z not in rs: fehler += 1
        nk = nullstelle(p_mul([p - n1, 0, 1], [p - n2, 0, 1], p), p, rng); faelle += 1
        if nk is not None: fehler += 1
    ok = fehler == 0
    gruen['K4'] = ok
    sag(f'K4 Nullstelle: {faelle} Faelle (p in (50, 500) und 2^127−1, je vier PK + eine NK), {fehler} Fehler ⇒ {"✅" if ok else "❌"}')

    # K5: `H23` = coefficients of H_23 (lowest first); `za`, `zb`, `zc` = number of cases (a), (b), (c); `fa`,
    #     `fb`, `fc` = failures in them
    H23 = [23375 ** 3, -5151296875, 3491750, 1]
    za, zb, zc, fa, fb, fc = 0, 0, 0, 0, 0, 0
    for p in [p for p in PR if p > 200]:
        f = [c % p for c in H23]
        n_lin = grad(linearteil(f, p))   # number of distinct roots of H_23 modulo p
        if cornacchia4(p, 23) is not None:
            za += 1; fa += n_lin != 3
        elif jacobi(-23 % p, p) == 1:
            zb += 1; fb += n_lin != 0
        elif jacobi(-23 % p, p) == -1:
            zc += 1; fc += n_lin != 1
    gruen['K5a'] = fa == 0 and za > 0
    gruen['K5b'] = fb == 0 and zb > 0
    gruen['K5c'] = fc == 0 and zc > 0
    sag(f'K5 H₂₃ (Atkin–Morain S. 43), Primzahlen 200 < p < 3000: (a) 4p = A² + 23B² loesbar: {za} Faelle, {fa} ohne drei Nullstellen '
        f'{"✅" if gruen["K5a"] else "❌"} · (b) Rest, aber keine Darstellung: {zb} Faelle, {fb} mit Nullstelle {"✅" if gruen["K5b"] else "❌"} · '
        f'(c) Nichtrest: {zc} Faelle, {fc} nicht genau eine Nullstelle {"✅" if gruen["K5c"] else "❌"}')

    alle = all(gruen.values())
    sag(f'ERGEBNIS: {"✅ alle Kontrollen gruen" if alle else "❌ " + str([k for k, v in gruen.items() if not v])} · {time.time() - t0:.1f} s')
    ERG.mkdir(exist_ok=True)
    text = '\n'.join(aus) + '\n'
    (ERG / 'e1_grundbausteine_output.txt').write_text(text, encoding='utf-8', newline='\n')
    sag(f'md5 Ausgabe {hashlib.md5(text.encode("utf-8")).hexdigest()[:8]} · md5 Skript {hashlib.md5(pathlib.Path(__file__).read_bytes()).hexdigest()[:8]}')
    return 0 if alle else 1

if __name__ == '__main__':
    sys.exit(main())
