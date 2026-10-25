# Part of the data deposit of 'Consecutive Powerful Numbers and Pell Equations' (B. Y. Sayim, 2026).
# License: Apache-2.0, see LICENSE and NOTICE.
#
# e5_teil2.py
# -*- coding: utf-8 -*-
# Own ECPP prover, part II: builds our own primality certificates for the 19 ECPP numbers of part II (48 to 1687 digits).
# Reads:  ergebnisse/w165_zertifikate/w165_auftrag.json (list of numbers with fields `i`, `zahl`, `stellen`);
#         uses e3_ecpp_schritt and e5_teil1.
# Writes: ecpp/ergebnisse/e5_cert_t2_<i>_eigen.txt (one certificate per number) and a log line per number
#         in ecpp/ergebnisse/e5_teil2_log.txt.
# Usage:  python e5_teil2.py <max_digits> <seconds_per_number> [<min_digits>]
#         only numbers with at most <max_digits> digits (and at least <min_digits>) are processed, in ascending order.
# Control: every certificate produced is re-checked with `pruefe_zertifikat` of w164 (the same checker w165 uses), and its first
#         entry must be the input number; the log records whether w164 accepted or rejected the chain.
# Expectation: the 15 numbers up to 357 digits are proved within minutes (part I: 234 digits in 43 s); 621 to 693 digits take
#         hours or may fail; 1191 and 1687 digits are not expected to finish in a single run.
import sys, time, random, json, pathlib
sys.stdout.reconfigure(encoding='utf-8')
HIER = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HIER))   # `HIER` = folder of this script
import e3_ecpp_schritt as e3
import e5_teil1 as t1
# `AUFTRAG2` = job file of part II (`auftrag` = order, job): list of {`i`, `zahl`, `stellen`}
AUFTRAG2 = HIER.parent / 'ergebnisse' / 'w165_zertifikate' / 'w165_auftrag.json'
def main():
    # `grenze_st` = digit limit, `grenze_s` = time limit per number in seconds, `ab_st` = optional minimum digit count,
    # `auf` = jobs sorted by digits
    grenze_st = int(sys.argv[1]); grenze_s = int(sys.argv[2]); ab_st = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    auf = sorted(json.load(open(AUFTRAG2, encoding='utf-8')), key=lambda x: int(x['stellen']))
    log = t1.ERG / 'e5_teil2_log.txt'
    for x in auf:
        if int(x['stellen']) > grenze_st or int(x['stellen']) < ab_st: continue
        i = x['i']; zahl = int(x['zahl']); rng = random.Random(2000 + int(i))
        # skip numbers that already have a certificate (only numbers still without one are processed)
        if (t1.ERG / f'e5_cert_t2_{i}_eigen.txt').exists(): print(f'Nr. {i}: Zertifikat liegt schon, uebersprungen', flush=True); continue
        print(f'E5 Teil II Nr. {i}: {len(str(zahl))} Stellen, Grenze {grenze_s} s', flush=True)
        # `kette` = certificate chain (None if not found), `infos` = per-step info, `status` = outcome text,
        # `sek` = seconds used
        kette, infos, status, sek = t1.lauf(i, zahl, grenze_s, rng)
        ok = False
        if kette is not None:
            ok = e3.pruefe_zertifikat([tuple(z) for z in kette]) is True and kette[0][0] == zahl
            (t1.ERG / f'e5_cert_t2_{i}_eigen.txt').write_text(repr(kette), encoding='utf-8', newline='\n')
        # `zeile` = log line; `lauf` (in e5_teil1) = run: builds the chain for one number
        zeile = (f'{time.strftime("%Y-%m-%d %H:%M")} T2-Nr. {i} ({len(str(zahl))} Stellen): {status}, {len(kette) if kette else len(infos)} Schritte, {sek:.0f} s, '
                 f'w164: {"✅ angenommen" if ok else ("❌ abgelehnt" if kette is not None else "kein Zertifikat")}')
        print(zeile, flush=True)
        with open(log, 'a', encoding='utf-8', newline='\n') as f: f.write(zeile + '\n')
if __name__ == '__main__':
    main()
