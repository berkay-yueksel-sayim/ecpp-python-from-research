# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# e5_teil1.py
# -*- coding: utf-8 -*-
# Purpose: the six large probable primes of Part I (111, 112, 194, 234, 245 and 477 digits, listed in
#   w164_zertifikate/w164_auftrag.json) proved with OUR OWN ECPP prover: our code (modules e3_ecpp_schritt and e4_besser)
#   generates its own certificates. Whether they are valid is decided by the independent checker (skripte/w164_ecpp_pruefer, and
#   PARI/GP's primecertisvalid), not by this script.
# Usage: python e5_teil1.py <numbers, e.g. 1,2,3> <seconds per number>
#   Reads ../ergebnisse/w164_zertifikate/w164_auftrag.json. Writes ergebnisse/e5_cert_<i>_eigen.txt (the chain, as a Python
#   literal) and appends a line to ergebnisse/e5_teil1_log.txt (paths relative to this folder).
# Step: taken from e4_besser: best of K candidates, smooth decomposition via gcd with a primorial
#   (the plain e3.schritt was too slow).
# Expectations, stated in advance:
#   E1  Nos. 1 and 2 (111/112 digits) are managed by the own build in minutes (conjecture from 60 digits taking 2-7 s and scaling
#       with about the fourth power of the number of digits; [assumption]).
#   E2  Every generated chain is accepted by the checker (otherwise it is an error in the own build, not in the checker's
#       verdict); the number in line 1 of the chain is the number from the order.
#   E3  No. 6 (477 digits) is NOT certain within the time limit with a pure Python implementation; no result there is not
#       an error but a measurement of the limit.

import sys, time, random, json, pathlib, hashlib
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import e3_ecpp_schritt as e3
import e4_besser as e4
# `ERG` = output folder; `AUFTRAG` = the order: list of the numbers to be proved
ERG = HIER / 'ergebnisse'; ERG.mkdir(exist_ok=True)
AUFTRAG = HIER.parent / 'ergebnisse' / 'w164_zertifikate' / 'w164_auftrag.json'

def lauf(i, zahl, grenze, rng):
    """Chain down to q <= 2^64. If a step finds no candidate (all three D lists empty), the chain goes back ONE stage and
    chooses a different candidate there (the failed q goes on the forbidden list of that stage); at most 40 backtracks
    per number. Returns (chain, infos, status, seconds); the chain is None on failure. (`lauf` = run, `grenze` = time limit)"""
    # `stufen` = stages (step, q, info); `verbot` = forbidden q per number N; `rueck` = number of backtracks; N = current number
    t0 = time.time(); stufen = []; verbot = {}; rueck = 0; N = zahl
    # `dl_klein`, `dl_gross`, `dl_riesig` = lists of discriminants D (small, large, huge), the larger ones built only when needed
    # (`hmax` = bound on the class number h, `dmax` = bound on |D|)
    dl_klein = e3.d_liste(); dl_gross = None; dl_riesig = None
    while N > 2 ** 64:
        if time.time() - t0 > grenze: return None, [s[2] for s in stufen], 'Zeitgrenze', time.time() - t0
        vb = verbot.get(N, set())
        r = e4.schritt_best(N, rng, dl_klein, verboten=vb)
        if r is None:
            if dl_gross is None: dl_gross = e3.d_liste(hmax=14, dmax=40000)
            r = e4.schritt_best(N, rng, dl_gross, verboten=vb)
        if r is None:
            if dl_riesig is None: dl_riesig = e3.d_liste(hmax=24, dmax=200000)
            r = e4.schritt_best(N, rng, dl_riesig, verboten=vb)
        if r is None:
            if not stufen or rueck >= 40: return None, [s[2] for s in stufen], f'kein Schritt bei {len(str(N))} Stellen (Rueckschritte {rueck})', time.time() - t0
            st, q, info = stufen.pop(); rueck += 1
            verbot.setdefault(st[0], set()).add(q); N = st[0]
            print(f'   Nr. {i}: kein Schritt bei {len(str(q))} Stellen — zurueck auf {len(str(N))} Stellen (Rueckschritt {rueck})', flush=True)
            continue
        st, q, info = r
        stufen.append((st, q, info))
        print(f'   Nr. {i}: Schritt {len(stufen)}: {len(str(N))} → {len(str(q))} Stellen, D = {info["D"]} (h = {info["h"]}), {time.time() - t0:.0f} s', flush=True)
        N = q
    return [s[0] for s in stufen], [s[2] for s in stufen], 'ok' + (f' (Rueckschritte {rueck})' if rueck else ''), time.time() - t0

def main():
    nummern = [int(x) for x in sys.argv[1].split(',')]; grenze = int(sys.argv[2])
    auftrag = {x['i']: x for x in json.load(open(AUFTRAG, encoding='utf-8'))}
    log = ERG / 'e5_teil1_log.txt'
    for i in nummern:
        # `zahl` = the number to prove, `rng` = random generator seeded per number
        zahl = int(auftrag[i]['zahl']); rng = random.Random(1000 + i)
        print(f'E5 Nr. {i}: {len(str(zahl))} Stellen, Grenze {grenze} s', flush=True)
        kette, infos, status, sek = lauf(i, zahl, grenze, rng)
        ok = False
        if kette is not None:
            # `ok` = `e3.pruefe_zertifikat` accepts the chain and its first entry is the given number (`kette` = chain)
            ok = e3.pruefe_zertifikat([tuple(z) for z in kette]) is True and kette[0][0] == zahl
            (ERG / f'e5_cert_{i}_eigen.txt').write_text(repr(kette), encoding='utf-8', newline='\n')
        zeile = (f'{time.strftime("%Y-%m-%d %H:%M")} Nr. {i} ({len(str(zahl))} Stellen): {status}, {len(kette) if kette else len(infos)} Schritte, {sek:.0f} s, '
                 f'w164: {"✅ angenommen" if ok else ("❌ abgelehnt" if kette is not None else "kein Zertifikat")}')
        print(zeile, flush=True)
        with open(log, 'a', encoding='utf-8', newline='\n') as f: f.write(zeile + '\n')

if __name__ == '__main__':
    main()
