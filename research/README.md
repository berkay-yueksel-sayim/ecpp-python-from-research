# The research version, unchanged

This folder is a copy of the ECPP part of the data deposit of the paper series *Consecutive Powerful Numbers and Pell
Equations: Triples, Pairs, and Further Observations*. The paths are the same as below `deposit/` in the data repository
[berkay-yueksel-sayim/erdos-364-365-powerful-numbers-pell](https://github.com/berkay-yueksel-sayim/erdos-364-365-powerful-numbers-pell),
and every file except this README has the same Git blob hash there (commit `fca8f016f720` of 13 October 2026, checked
on 23 October 2026). Nothing here has been
edited since that publication. (Before it, the deposit's scripts had internal working notes removed and a few printed
messages reworded, in the scripts and their stored outputs; no computed value changed, see `deposit/DATA_README.md` in
the data repository.) The package in
`ecpp_research/` is the cleaned-up version, and its tests load the checker below unchanged and compare the two.

The scripts are commented in English; many identifiers are German (`pruefe_zertifikat` = check certificate, `schritt` =
step, `ergebnisse` = results, `Teil` = Part).

## What is here

| Path | What it is |
|---|---|
| `skripte/w164_ecpp_pruefer_2026-10-02.py` | The checker that verified every elliptic curve certificate in the papers. Run on its own, it checks the six certificates of Part I and the control certificate. |
| `ergebnisse/w164_zertifikate/` | Part I: the PARI/GP 2.17.2 certificates `cert_1.txt` … `cert_6.txt` (111 to 477 digits), the positive control `cert_pk.txt` (PARI/GP's certificate for nextprime(10^40), read by w164), the numbers (`w164_auftrag.json`), the PARI job (`w164_ecpp.gp`) and its run time (`zeit.txt`). |
| `ergebnisse/w165_zertifikate/` | Part II: the PARI/GP certificates `cert_1.txt` … `cert_19.txt` (48 to 1687 digits), the numbers, the PARI jobs and run times. |
| `ergebnisse/w164_ecpp_pruefer_*`, `ergebnisse/w165_teil2_primbeweise_*` | Output and result file of the two runs. |
| `ecpp/e1_grundbausteine.py` … `e5_teil2.py` | Our own prover, written as a second route: building blocks (e1), Hilbert class polynomials on two routes (e2), one ECPP step (e3), best of several candidates (e4, with `e4_kontrolle.py`), the chains for the numbers of Part I and Part II (e5). |
| `ecpp/pari_pruefen.gp` | The third route: PARI's built-in verifier `primecertisvalid` applied to our own certificates. |
| `ecpp/ergebnisse/` | The outputs of e1, e2 and e3 (the chain part of e3 is labelled E4 in its output; `e4_besser.py` and `e4_kontrolle.py` write no file), the run logs of e5, our own certificates (`e3_cert_pk_1e40.txt`, `e5_cert_<i>_eigen.txt` for Part I, `e5_cert_t2_<i>_eigen.txt` for Part II), and `e5_pari_pruefung_output.txt` (PARI accepts 23 of 23). |

**Not included:** `w165_teil2_primbeweise_2026-10-02.py`, the script that ran PARI/GP for Part II and wrote
`w165_zertifikate/`. It imports another script of the deposit (w144, which rebuilds the numbers) and calls PARI/GP in a
container, so it does not run on its own; it is in the data repository under `deposit/skripte/`. Its checking functions are
verbatim copies of those in w164.

## A gap in w164, found after publication

w164 is not a sound checker for composite numbers. Its group arithmetic decides each case (is a point the point at
infinity, are two points equal) modulo N only. When N is composite, a step can look fine modulo N while the points are
degenerate modulo a prime divisor of N, and w164 then accepts certificates for some composite numbers; three such
certificates are in `../tests/test_soundness.py` (for example one for 10196873 = 1637 · 6229). The corrected checker in
`../ecpp_research/` refuses every such step (see the main README, condition 4).

This does not change any result of the papers: all 25 PARI/GP certificates and all 23 of our own are valid, which the
corrected checker confirms, and PARI's `primecertisvalid` accepts all 48 as well. But w164 should not be used to check
new certificates; use the package.

## The routes

| | checked by w164 | checked by the corrected checker (`ecpp_research`) | checked by PARI's `primecertisvalid` |
|---|---|---|---|
| 25 PARI/GP certificates (6 for Part I, 19 for Part II) | 25 of 25 | 25 of 25 | 25 of 25 (`../checks/pari_all_25_output.txt`) |
| 23 certificates from our own prover (6 for Part I, 17 for Part II) | 23 of 23 | 23 of 23 | 23 of 23 (`ecpp/ergebnisse/e5_pari_pruefung_output.txt`) |

The two control certificates for nextprime(10^40) (`cert_pk.txt` from PARI/GP, `e3_cert_pk_1e40.txt` from our prover)
are accepted by w164 and the corrected checker; PARI checked `cert_pk.txt` as the positive control of
`../checks/pari_all_25.gp`. The PARI checks of the 25 PARI certificates and the checks by the corrected checker
were run for this repository (October 2026), after the papers were published.

Our own prover produced no certificate for the primes of `w165_zertifikate/cert_16.txt` (1687 digits: no first step found
within 344 seconds, `ecpp/ergebnisse/e5_teil2_log.txt`) and `cert_6.txt` (1191 digits: the run left no result); these two
rest on the PARI certificates alone. (The certificate files are numbered by computation job; the papers do not number
these primes.)

## Running it

```bash
cd research/skripte
python w164_ecpp_pruefer_2026-10-02.py
```

The checker needs only the Python standard library and calls no other program. The scripts in `ecpp/` run the same way
from their own folder (`python e1_grundbausteine.py`, …); e2 to e5 need `mpmath`. Both w164 and the scripts in `ecpp/`
write into their `ergebnisse/` folder (w164 rewrites `w164_ecpp_pruefer_result.json` and `w164_ecpp_pruefer_output.txt`),
so run them in a copy if you want to keep the files as published. `pari_pruefen.gp` is a script for PARI/GP 2.17.2; its
header says where it expects the certificate files.

The tests of the package load the functions of w164 from its source text, up to the line where its main run starts
(`../tests/conftest.py`), and compare them with the package. They rely on w164 staying unchanged.

## License

The Python scripts and `ecpp/pari_pruefen.gp` are licensed under the Apache License, Version 2.0 (`../LICENSE`,
`../NOTICE`); every file in `ergebnisse/` and `ecpp/ergebnisse/`, including the PARI/GP job files `w164_ecpp.gp`,
`w165_ecpp.gp` and `w165_rest.gp`, is licensed under CC BY 4.0 (`../LICENSE-data`), as in the data repository.
