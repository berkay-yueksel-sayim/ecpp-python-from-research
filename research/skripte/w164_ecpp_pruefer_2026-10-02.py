# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# w164_ecpp_pruefer_2026-10-02.py
# -*- coding: utf-8 -*-
# Our own checker for the ECPP certificates (PARI primecert format) of the six large probable primes.
#
# Reads:    ../ergebnisse/w164_zertifikate/cert_pk.txt (positive control), cert_<i>.txt (the six certificates, one per job) and
#           w164_auftrag.json (the job list: m, d, digit count, number, certificate index i).
# Writes:   ../ergebnisse/w164_ecpp_pruefer_result.json and ../ergebnisse/w164_ecpp_pruefer_output.txt
# Usage:    python w164_ecpp_pruefer_2026-10-02.py   (no arguments; pure Python, no outside library)
# Controls: E1 (positive), E2 (negative), E3 (all six certificates accepted), see below.
#
# BACKGROUND: after the Pocklington step (w163), six numbers of 111, 112, 194, 234, 245 and 477 digits (inventory w162)
#   remained. PARI/GP 2.17.2 (run in a container) generated their certificates (w164_zertifikate/cert_i.txt). The generating
#   program is NOT trusted: every certificate is checked by this code. A second, independent check is done by PARI's own
#   verifier (pari_pruefen.gp).
# THEOREM USED (a theorem, not code): Goldwasser-Kilian 1986 / Atkin-Morain 1993. Let N > 1 with gcd(N, 6) = 1, let
#   E: y² = x³ + a x + b with gcd(4a³ + 27b², N) = 1, m = s·q with q prime and q > (N^(1/4) + 1)², and let P be a point with
#   m·P = O and s·P ≠ O modulo every prime divisor p of N. Then N is prime. (For p | N, Q = s·P has order q modulo p, so
#   q ≤ #E(F_p) ≤ (√p + 1)²; q > (N^(1/4) + 1)² then gives p > √N.) Certificate format as described in the PARI documentation
#   (documentation only, no source code): C[i] = [N_i, t_i, s_i, a_i, P_i], m_i = N_i + 1 − t_i, q_i = m_i/s_i = N_{i+1}, the
#   last q_l ≤ 2^64 prime; b_i is recovered from the point: b_i = y_i² − x_i³ − a_i x_i.
# COMPUTATION: projective coordinates (X : Y : Z) modulo N, without inverses. The case distinctions (equal/opposite points) are
#   made modulo N; if one holds modulo N, it holds modulo every p. If it holds only modulo p, the formula degenerates modulo p to
#   (0 : 0 : 0) or to O; both make the Z coordinate zero modulo p and fail the checks below. Only the following is assumed:
#   (i)  Q = s·P with gcd(Z_Q, N) = 1  (Q is a genuine point modulo every p), and
#   (ii) (q − 1)·Q = −Q with gcd(Z, N) = 1 and the same projective class (then q·Q = O modulo every p, without degeneration).
#   The last q_l ≤ 2^64 < psi_12: deterministic test with the first 12 prime bases.
# EXPECTATIONS (set before the run):
#   E1  PK: PARI's certificate for nextprime(10^40) (read from cert_pk.txt) is accepted.
#   E2  NK: (a) N replaced by N + 2, (b) t shifted by 1, (c) a point coordinate corrupted, (d) s doubled; each is rejected.
#   E3  All six certificates are accepted; the number in line 1 of each certificate is exactly the number from w162
#       (job file w164_auftrag.json).
import sys, json, time, ast, pathlib
from math import gcd, isqrt
sys.stdout.reconfigure(encoding='utf-8')
# `HIER` = this directory, `ERG` = results, `Z` = certificates
HIER = pathlib.Path(__file__).resolve().parent; ERG = HIER.parent/'ergebnisse'; Z = ERG/'w164_zertifikate'
aus = []   # `aus` = lines of the output file
def sag(s=''):   # `sag` = print a line and record it for the output file
    print(s, flush=True); aus.append(s)
BASEN = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)   # `BASEN` = the 12 prime bases of the deterministic Miller-Rabin test
# `mr12` = deterministic Miller-Rabin test with the 12 bases above (correct for n < psi_12, in particular n <= 2^64)
def mr12(n):
    if n < 2: return False
    for p in BASEN:
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in BASEN:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

O = None                                      # the point at infinity (also any point with Z ≡ 0 modulo N, see `ist_O`)
def ist_O(P, N): return P is None or P[2] % N == 0   # `ist_O` = test for the point at infinity
def addiere(P, Q, a, N):   # `addiere` = addition P + Q on y² = x³ + a x + b in projective coordinates modulo N
    if ist_O(P, N): return Q
    if ist_O(Q, N): return P
    X1, Y1, Z1 = P; X2, Y2, Z2 = Q
    u = (Y2 * Z1 - Y1 * Z2) % N; v = (X2 * Z1 - X1 * Z2) % N
    if v == 0:
        return verdopple(P, a, N) if u == 0 else O
    w = Z1 * Z2 % N; v2 = v * v % N; v3 = v2 * v % N
    A = (u * u * w - v3 - 2 * v2 * X1 * Z2) % N
    return (v * A % N, (u * (v2 * X1 * Z2 - A) - v3 * Y1 * Z2) % N, v3 * w % N)
def verdopple(P, a, N):   # `verdopple` = doubling 2P in projective coordinates modulo N
    if ist_O(P, N): return O
    X, Y, Zk = P
    if Y % N == 0: return O
    w = (a * Zk * Zk + 3 * X * X) % N; s = Y * Zk % N; B = X * Y * s % N; h = (w * w - 8 * B) % N
    return (2 * h * s % N, (w * (4 * B - h) - 8 * Y * Y * s * s) % N, 8 * s * s * s % N)
def mal(k, P, a, N):   # `mal` = scalar multiplication k·P by double-and-add over the binary digits of k
    R = O
    for bit in bin(k)[2:]:
        R = verdopple(R, a, N)
        if bit == '1': R = addiere(R, P, a, N)
    return R

def pruefe_schritt(Ni, t, s, a, P):
    # `pruefe_schritt` = check one certificate step (N_i, t, s, a, P); returns q = m/s (an int) if the step is valid,
    # otherwise a string naming the failed condition
    x, y = P
    if Ni < 2 or gcd(Ni, 6) != 1: return 'ggT(N, 6) ≠ 1'
    if t * t >= 4 * Ni: return 't² ≥ 4N'
    m = Ni + 1 - t
    if s <= 0 or m % s: return 's teilt m nicht'
    q = m // s
    r = isqrt(isqrt(Ni)) + 1                  # ≥ N^(1/4)
    if q <= (r + 1) ** 2: return 'q zu klein'   # requires q > (⌊N^(1/4)⌋ + 2)², stricter than necessary
    b = (y * y - x * x * x - a * x) % Ni
    if gcd((4 * a ** 3 + 27 * b * b) % Ni, Ni) != 1: return 'Kurve singulaer modulo einem Teiler'
    Pp = (x % Ni, y % Ni, 1)
    Q = mal(s, Pp, a, Ni)
    if Q is None or gcd(Q[2], Ni) != 1: return 's·P ist modulo einem Teiler O'
    R = mal(q - 1, Q, a, Ni)
    if R is None or gcd(R[2], Ni) != 1: return '(q−1)·Q entartet'
    X1, Y1, Z1 = R; X2, Y2, Z2 = Q
    if (X1 * Z2 - X2 * Z1) % Ni or (Y1 * Z2 + Y2 * Z1) % Ni: return '(q−1)·Q ≠ −Q'
    return q
def pruefe_zertifikat(C):
    # `pruefe_zertifikat` = check a whole certificate: every step valid, N_{i+1} = q_i, last q small and prime;
    #   returns True or an error string
    erwartet = None   # `erwartet` = the value N_{i+1} that the next step must have (q of the previous step)
    for Ni, t, s, a, P in C:
        if erwartet is not None and Ni != erwartet: return 'Kette bricht (N_{i+1} ≠ q_i)'
        q = pruefe_schritt(Ni, t, s, a, P)
        if isinstance(q, str): return q
        erwartet = q
    if erwartet is None or erwartet > 2 ** 64 or not mr12(erwartet): return 'letztes q nicht klein und prim'
    return True

def lies(p):   # `lies` = read a certificate file (a Python literal: list of steps) and return the steps as tuples
    return [tuple(z) for z in ast.literal_eval(p.read_text(encoding='utf-8').strip())]

# E1 PK: PARI's certificate for nextprime(10^40), read from cert_pk.txt; if the file is missing, the PK is reported as skipped
#   and the run exits with status 1
pk = Z / 'cert_pk.txt'
if pk.exists():
    e1 = pruefe_zertifikat(lies(pk)) is True
    sag(f'PK E1 {"✅" if e1 else "❌"}  PARI-Zertifikat fuer nextprime(10^40) angenommen'); assert e1
else:
    sag('PK E1 ⚠️  cert_pk.txt fehlt — PK uebersprungen'); raise SystemExit(1)
# E2 NK
C0 = lies(pk)
def kopie(): return [list(z) for z in C0]   # `kopie` = fresh mutable copy of the control certificate
n1 = kopie(); n1[0][0] += 2                                  # (a) N replaced by N + 2
n2 = kopie(); n2[0][1] += 1                                  # (b) t shifted by 1
n3 = kopie(); n3[0][4] = [n3[0][4][0], n3[0][4][1] + 1]      # (c) y coordinate of the point corrupted
n4 = kopie(); n4[0][2] *= 2                                  # (d) s doubled
nk = [pruefe_zertifikat([tuple(z) for z in n]) for n in (n1, n2, n3, n4)]   # `nk` = results for the four corrupted copies
e2 = all(r is not True for r in nk)
sag(f'NK E2 {"✅" if e2 else "❌"}  verfaelscht abgelehnt: N+2 → „{nk[0]}" · t+1 → „{nk[1]}" · Punkt → „{nk[2]}" · 2s → „{nk[3]}"'); assert e2

auftrag = json.load(open(Z / 'w164_auftrag.json', encoding='utf-8'))   # `auftrag` = job list
ergebnisse = []   # `ergebnisse` = one result record per job
for x in auftrag:
    C = lies(Z / f'cert_{x["i"]}.txt'); t0 = time.time()
    ok_zahl = C[0][0] == int(x['zahl'])   # the certificate must belong to the number of the job
    r = pruefe_zertifikat(C) if ok_zahl else 'Zertifikat gehoert zu einer anderen Zahl'
    ok = r is True
    ergebnisse.append(dict(m=x['m'], d=x['d'], stellen=x['stellen'], zahl=x['zahl'], bewiesen=ok, schritte=len(C), pruefung=('ok' if ok else r),
                           sek=round(time.time() - t0, 1)))
    sag(f'  ({x["m"]}, {x["d"]}) {x["stellen"]:>3} St.: ' + (f'✅ bewiesen prim (ECPP, {len(C)} Schritte), eigener Pruefer' if ok else f'🔴 abgelehnt: {r}')
        + f'  {time.time() - t0:.1f} s')
n_ok = sum(1 for e in ergebnisse if e['bewiesen'])   # `n_ok` = number of accepted certificates
sag(f'E3 {"✅" if n_ok == len(ergebnisse) else "❌"}  {n_ok} von {len(ergebnisse)} ECPP-Zertifikaten angenommen')
res = dict(skript=pathlib.Path(__file__).name, datum=time.strftime('%Y-%m-%d %H:%M'), werkzeug='PARI/GP 2.17.2 (Debian pari-gp 2.17.2-1, Image dkr-pari:1)',
           n=len(ergebnisse), bewiesen=n_ok, ergebnisse=ergebnisse)
(ERG/'w164_ecpp_pruefer_result.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
(ERG/'w164_ecpp_pruefer_output.txt').write_text(f'w164 · {time.strftime("%Y-%m-%d %H:%M")}\n' + '\n'.join(aus) + '\n', encoding='utf-8')
sag('Ergebnis: w164_ecpp_pruefer_result.json')
