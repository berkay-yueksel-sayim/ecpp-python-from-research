# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# e3_ecpp_schritt.py
# -*- coding: utf-8 -*-
# Own ECPP prover, E3/E4: one ECPP step from our own hands, and the chain down to q ≤ 2^64.
#
# Reads:    e1_grundbausteine.py and e2_klassenpolynome.py (same folder, imported), and ../skripte/w164_ecpp_pruefer_2026-10-02.py
#           (its functions are loaded from the source text, without its main run).
# Writes:   ergebnisse/e3_cert_pk_1e40.txt (our certificate chain for nextprime(10^40)) and ergebnisse/e3_ecpp_schritt_output.txt
#           (both relative to this file), and printed output.
# Usage:    python e3_ecpp_schritt.py   (no arguments; pure Python, no outside library); exit code 0 if all controls are green.
# Controls: K1 to K5 (listed below); PK = positive control, NK = negative control.
#
# BASIS (theorems, no pseudocode): Atkin–Morain 1993 — Theorem 5.2 / Cor. 5.1 (the checking criterion; the checker w164
#   implements it), §5.3 (procedure: D with N = ππ' in Q(√−D), candidates m = N + 1 ∓ A, m = s·q with q > (N^{1/4} + 1)²,
#   curve from a root of H_D, point, check), §8.6.2 Eq. (39) (the curve y² = x³ + 3kc²x + 2kc³ with k = j/(1728 − j) has
#   invariant j; c switches between the curve and its twist), §8.6.3 (point without a square root: λ = x₀³ + a x₀ + b, P =
#   (λx₀, λ²) on Y² = X³ + aλ²X + bλ³). Building blocks from e1/e2. The cases D = 3, 4 (j = 0, 1728, six resp. four
#   candidates m) are LEFT OUT here: we choose D ≥ 7.
# DECIDER: `w164` alone (functions `pruefe_schritt` / `pruefe_zertifikat`, loaded from the source text of the deposit,
#   without its main run). This script only SEARCHES; what it finds counts only once w164 accepts it. Therefore an error
#   here is never a false proof.
# EXPECTATIONS (set before the run):
#   K1  PK: nextprime(10^40): our chain is accepted by w164 (the same number that PARI's cert_pk.txt
#       proves; the chains may differ).
#   K2  PK: three random (probable) primes each with 30, 45 and 60 digits: chains accepted, within the time limit.
#   K3  NK: composite N (products of two 20-digit primes): NO accepted step arises.
#   K4  NK: an accepted step in which t is shifted by 1, s is doubled, or the point is corrupted fails at w164.
#   K5  Reference values Atkin–Morain pp. 59/60: for p = 17401, D = 8 and p = 107, D = 7 the step finds a curve with m = p + 1 −
#       t, and the point counts 17200 and 88 given there occur as one of the two candidates m.
import sys, math, time, random, pathlib, hashlib, json
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent   # `HIER` = this directory
sys.path.insert(0, str(HIER))
from e1_grundbausteine import jacobi, cornacchia4, nullstelle, p_ggt, linearteil, grad
from e2_klassenpolynome import klassenpolynom, reduzierte_formen, fundamental
ERG = HIER / 'ergebnisse'   # `ERG` = results directory
W164 = HIER.parent / 'skripte' / 'w164_ecpp_pruefer_2026-10-02.py'   # `W164` = path of the checker script

# ---------- the checker: load the w164 functions without its main run (which starts at the line `pk = Z /
#   'cert_pk.txt'`) ----------
def lade_w164():   # `lade_w164` = load w164: executes the source text up to the main run and returns the namespace
    src = W164.read_text(encoding='utf-8')
    kopf = src.split("pk = Z / 'cert_pk.txt'")[0]
    ns = {'__file__': str(W164), '__name__': 'w164_funktionen'}
    exec(compile(kopf, str(W164), 'exec'), ns)
    return ns
w164 = lade_w164()
pruefe_schritt, pruefe_zertifikat, mal, ist_O = w164['pruefe_schritt'], w164['pruefe_zertifikat'], w164['mal'], w164['ist_O']

# ---------- helpers ----------
def primzahlen_bis(n):   # sieve of Eratosthenes: all primes <= n
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, math.isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
KLEINE = primzahlen_bis(200000)   # `KLEINE` = small primes used for trial division
# strong probable-prime test with random bases — only for CANDIDATE SELECTION; the chain supplies the proof
def wprim(n, rng, runden=24):
    if n < 2: return False
    for p in KLEINE[:50]:
        if n % p == 0: return n == p
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for _ in range(runden):
        a = rng.randrange(2, n - 1); x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
# Pollard rho (Brent variant) — own standard implementation; returns a factor or None
def rho_brent(n, rng, max_iter=200000):
    if n % 2 == 0: return 2
    y, c, m = rng.randrange(1, n), rng.randrange(1, n), 128
    g, r, q = 1, 1, 1; it = 0
    while g == 1:
        x = y
        for _ in range(r): y = (y * y + c) % n
        k = 0
        while k < r and g == 1:
            ys = y
            for _ in range(min(m, r - k)):
                y = (y * y + c) % n; q = q * abs(x - y) % n
            g = math.gcd(q, n); k += m; it += m
            if it > max_iter: return None
        r *= 2
    if g == n:
        while True:
            ys = (ys * ys + c) % n; g = math.gcd(abs(x - ys), n)
            if g > 1: break
    return None if g == n else g
def zerlege_m(m, N, rng, rho_iter):   # `zerlege_m` = split the order m into s and a probable prime q (see the docstring)
    """m = s·q with q probably prime and q > (N^(1/4)+1)² (the bound of w164: q > (⌊N^(1/4)⌋ + 2)², stricter). Otherwise None."""
    schranke = (math.isqrt(math.isqrt(N)) + 2) ** 2   # `schranke` = lower bound for q
    s, rest = 1, m
    for p in KLEINE:
        if p * p > rest: break
        while rest % p == 0: rest //= p; s *= p
    while rest > 1:
        if rest <= schranke: return None                   # the large prime part has become too small
        if wprim(rest, rng): return (s, rest) if s * rest == m else None
        f = rho_brent(rest, rng, rho_iter)
        if f is None or f in (1, rest): return None
        klein, gross = min(f, rest // f), max(f, rest // f)
        # the smaller factor (even if composite) goes into s; continue with the large one
        s *= klein; rest = gross
    return None

# ---------- discriminants and class polynomials ----------
# `d_liste` = fundamental discriminants D >= 7 below dmax with class number h <= hmax, sorted by (h, D)
def d_liste(hmax=8, dmax=6000):
    L = []
    for D in range(7, dmax):
        if D % 4 not in (0, 3) or not fundamental(D): continue
        h = len(reduzierte_formen(D, abbruch=hmax))
        if h <= hmax: L.append((h, D))
    L.sort()
    return [D for h, D in L]
HD = {}   # `HD` = cache of the class polynomials H_D
def hd(D):   # class polynomial H_D (first return value of `klassenpolynom`), cached
    if D not in HD: HD[D] = klassenpolynom(D)[0]
    return HD[D]

# ---------- one step ----------
def schritt(N, rng, DL, rho_iter=200000, versuche_punkt=12):
    """Searches for N (probably prime, N > 2^64) a step [N, t, s, a, [x, y]] that w164 accepts. Returns (step, q, info) or None."""
    for D in DL:   # `DL` = list of discriminants D to try
        if jacobi(-D % N, N) != 1: continue
        darst = cornacchia4(N, D)   # `darst` = representation 4N = A² + D·B²
        if darst is None: continue
        A, B = darst
        kandidaten = [N + 1 - A, N + 1 + A]   # `kandidaten` = candidate group orders m
        zerl = {}   # `zerl` = candidate m -> (s, q) for those m that split suitably
        for m in kandidaten:
            z = zerlege_m(m, N, rng, rho_iter)
            if z: zerl[m] = z
        if not zerl: continue
        f = [c % N for c in hd(D)]
        j = nullstelle(f, N, rng)   # `j` = a root of H_D modulo N (a j-invariant)
        if j is None or j == 0 or j == 1728 % N: continue
        k = j * pow((1728 - j) % N, -1, N) % N
        a0, b0 = 3 * k % N, 2 * k % N   # curve y² = x³ + a0 x + b0 with invariant j (c = 1)
        for _ in range(versuche_punkt):
            x0 = rng.randrange(1, N)
            lam = (x0 * x0 * x0 + a0 * x0 + b0) % N   # `lam` = λ, makes a point without a square root possible
            if lam == 0: continue
            a, b = a0 * lam * lam % N, b0 * lam * lam * lam % N
            P = (lam * x0 % N, lam * lam % N)
            for m, (s, q) in zerl.items():
                t = N + 1 - m
                r = pruefe_schritt(N, t, s, a, P)
                if isinstance(r, int) and r == q:
                    return [N, t, s, a, [P[0], P[1]]], q, dict(D=D, A=A, B=B, m=m, h=len(reduzierte_formen(D)))
    return None

# `zertifikat` = builds the chain of steps from N down to q <= 2^64; returns (chain, infos, status)
def zertifikat(N, rng, DL, zeitgrenze=600, **kw):
    kette, infos, t0 = [], [], time.time()
    while N > 2 ** 64:
        if time.time() - t0 > zeitgrenze: return None, infos, 'Zeitgrenze'
        r = schritt(N, rng, DL, **kw)
        if r is None: return None, infos, f'kein Schritt fuer {len(str(N))}-stellige Zahl'
        st, q, info = r
        kette.append(st); infos.append(info); N = q
    return kette, infos, 'ok'

# ================= Controls =================
aus = []   # `aus` = lines of the output file
def sag(s=''):   # `sag` = print a line and record it for the output file
    print(s, flush=True); aus.append(s)
def naechste_prim(n, rng):                         # smallest probable prime > n
    # make n odd and larger (an earlier form `n += 1 - (n % 2 == 0)` left an even n even: endless loop)
    n = n + 1 if n % 2 == 0 else n + 2
    while not wprim(n, rng): n += 2
    return n

def main():
    # `gruen` = status ("green" = passed) of each control; `DL` = the discriminants to try
    t0 = time.time(); rng = random.Random(20261003); gruen = {}
    DL = d_liste()
    sag(f'E3/E4 — eigener ECPP-Schritt und eigene Kette · Pruefer: w164-Funktionen aus {W164.name} · {len(DL)} Diskriminanten mit h ≤ 8 (D ≥ 7)')

    # K5 first (cheap): reference values from Atkin–Morain
    k5 = []
    for p, D, soll in ((17401, 8, 17200), (107, 7, 88)):   # `soll` = the point count given in the paper
        darst = cornacchia4(p, D); kand = None if darst is None else (p + 1 - darst[0], p + 1 + darst[0])
        k5.append(kand is not None and soll in kand)
        sag(f'K5 p = {p}, D = {D}: 4p = A² + {D}B² → {darst}; Kandidaten m = {kand}; {soll} dabei: {"✅" if k5[-1] else "❌"}')
    gruen['K5'] = all(k5)

    # K1: nextprime(10^40) — chain, then w164
    N = naechste_prim(10 ** 40, rng)
    kette, infos, status = zertifikat(N, rng, DL, zeitgrenze=600)   # `kette` = our chain of steps, `infos` = data of each step
    ok1 = kette is not None and pruefe_zertifikat([tuple(z) for z in kette]) is True
    gruen['K1'] = ok1
    sag(f'K1 nextprime(10^40) = {N}: {status}, {len(kette) if kette else 0} Schritte, D = {[i["D"] for i in infos]}, '
        f'w164: {"✅ angenommen" if ok1 else "❌"} · {time.time() - t0:.1f} s')
    if kette:
        ERG.mkdir(exist_ok=True)
        (ERG / 'e3_cert_pk_1e40.txt').write_text(repr(kette), encoding='utf-8', newline='\n')

    # K4: falsifications of the first step (t + 1, s doubled, point corrupted); each must be rejected by w164
    if kette:
        st = kette[0]
        v = [list(st) for _ in range(3)]
        v[0][1] += 1; v[1][2] *= 2; v[2][4] = [st[4][0], st[4][1] + 1]
        r = [pruefe_schritt(*x[:4], x[4]) for x in v]
        gruen['K4'] = all(isinstance(x, str) for x in r)
        sag(f'K4 verfaelscht (t+1, 2s, Punkt): {r} ⇒ {"✅ alle abgelehnt" if gruen["K4"] else "❌"}')
    else:
        gruen['K4'] = False

    # K2: random primes with 30/45/60 digits, three each
    k2, zeiten = [], []   # `k2` = pass/fail per prime; `zeiten` = run time per prime in seconds
    for stellen in (30, 45, 60):
        for _ in range(3):
            N = naechste_prim(rng.randrange(10 ** (stellen - 1), 10 ** stellen), rng)
            t1 = time.time(); kette, infos, status = zertifikat(N, rng, DL, zeitgrenze=600)
            ok = kette is not None and pruefe_zertifikat([tuple(z) for z in kette]) is True
            k2.append(ok); zeiten.append(round(time.time() - t1, 1))
            sag(f'K2 {stellen} St.: {status}, {len(kette) if kette else 0} Schritte, w164 {"✅" if ok else "❌"}, {zeiten[-1]} s')
    gruen['K2'] = all(k2)

    # K3: composite numbers — no accepted step
    p1, p2 = naechste_prim(rng.randrange(10 ** 19, 10 ** 20), rng), naechste_prim(rng.randrange(10 ** 19, 10 ** 20), rng)
    carmichael = 1729  # unused: too small (≤ 2^64) for a chain; instead: 3 composite numbers of 40 digits
    zs = [p1 * p2] + [naechste_prim(rng.randrange(10 ** 19, 10 ** 20), rng) * naechste_prim(rng.randrange(10 ** 19, 10 ** 20), rng) for _ in range(2)]   # `zs` = the composite test numbers
    k3 = []   # `k3` = for each composite number: True if no step was found
    for Nz in zs:
        t1 = time.time()
        r = schritt(Nz, rng, DL[:12], rho_iter=50000, versuche_punkt=4)
        k3.append(r is None)
        sag(f'K3 zusammengesetzt ({len(str(Nz))} St.): {"✅ kein angenommener Schritt" if r is None else "❌ Schritt angenommen"} · {time.time() - t1:.1f} s')
    gruen['K3'] = all(k3)

    alle = all(gruen.values())
    sag(f'ERGEBNIS: {"✅ alle Kontrollen gruen" if alle else "❌ " + str([k for k, v in gruen.items() if not v])} · {time.time() - t0:.1f} s')
    ERG.mkdir(exist_ok=True)
    text = '\n'.join(aus) + '\n'
    (ERG / 'e3_ecpp_schritt_output.txt').write_text(text, encoding='utf-8', newline='\n')
    sag(f'md5 Ausgabe {hashlib.md5(text.encode("utf-8")).hexdigest()[:8]} · md5 Skript {hashlib.md5(pathlib.Path(__file__).read_bytes()).hexdigest()[:8]}')
    return 0 if alle else 1

if __name__ == '__main__':
    sys.exit(main())
