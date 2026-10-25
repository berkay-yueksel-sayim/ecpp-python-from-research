# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# e4_kontrolle.py
# Control run for the ECPP prover variant e4_besser.py (expectations F1-F3 below, taken from the header of e4_besser.py).
# Reads:  ergebnisse/w164_zertifikate/w164_auftrag.json (the list of numbers to be certified; entry i = 2 is used).
# Writes: nothing, prints one line per run and a final summary line.
# Usage:  python e4_kontrolle.py   (no arguments; fixed seed for the random numbers)
# Controls: every certificate chain built by e4_besser is checked by the independent checker w164 (via e3.pruefe_zertifikat).
#   F1: chains for nextprime(10^40) and for random primes of 30, 45, 60 digits are accepted by the checker.
#   F2: entry 2 of the order file (112 digits) is certified in under 60 s and the first step starts at that number.
#   F3: the median decrease in digit count per step for that number is at least 4.
import sys, time, random, json, pathlib, statistics
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HIER))
import e3_ecpp_schritt as e3, e4_besser as e4
rng = random.Random(20261003); DL = e3.d_liste()   # `DL` = list of candidate discriminants from e3
def kette(N, grenze=300):   # `kette` = chain: repeat the "best-of-K" step on N until N <= 2^64; `grenze` = time limit in s
    t0 = time.time(); k = []; red = []   # `k` = steps of the certificate; `red` = digits removed per step
    while N > 2 ** 64:
        if time.time() - t0 > grenze: return None, red
        r = e4.schritt_best(N, rng, DL)
        if r is None: return None, red
        st, q, info = r; k.append(st); red.append(len(str(N)) - len(str(q))); N = q
    return k, red
t0 = time.time(); alle = []   # `alle` = list of acceptance flags of the F1 runs
N = e3.naechste_prim(10 ** 40, rng)
k, red = kette(N); ok = k is not None and e3.pruefe_zertifikat([tuple(z) for z in k]) is True; alle.append(ok)
print(f'F1 nextprime(10^40): {len(k) if k else 0} Schritte, w164 {"✅" if ok else "❌"}, {time.time() - t0:.1f} s, Abnahme je Schritt {red}', flush=True)
for stellen in (30, 45, 60):
    for _ in range(3):
        N = e3.naechste_prim(rng.randrange(10 ** (stellen - 1), 10 ** stellen), rng); t1 = time.time()
        k, red = kette(N); ok = k is not None and e3.pruefe_zertifikat([tuple(z) for z in k]) is True; alle.append(ok)
        print(f'F1 {stellen} St.: {len(k) if k else 0} Schritte, w164 {"✅" if ok else "❌"}, {time.time() - t1:.1f} s', flush=True)
a = json.load(open(HIER.parent / 'ergebnisse/w164_zertifikate/w164_auftrag.json', encoding='utf-8'))
N = int([x for x in a if x['i'] == 2][0]['zahl']); t1 = time.time()
k, red = kette(N, 120); ok = k is not None and e3.pruefe_zertifikat([tuple(z) for z in k]) is True and k[0][0] == N
dt = time.time() - t1
print(f'F2 Nr. 2 (112 Stellen): {len(k) if k else 0} Schritte in {dt:.1f} s, w164 {"✅" if ok else "❌"} · F2 {"✅" if ok and dt < 60 else "❌"} · '
      f'F3 Median Stellenabnahme je Schritt: {statistics.median(red) if red else "—"} {"✅" if red and statistics.median(red) >= 4 else "❌"}', flush=True)
print('ERGEBNIS:', '✅ alle' if all(alle) and ok and dt < 60 else '❌ nicht alle', f'· {time.time() - t0:.0f} s')
