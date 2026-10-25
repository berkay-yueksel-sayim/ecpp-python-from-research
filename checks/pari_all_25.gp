\\ SPDX-License-Identifier: Apache-2.0
/* pari_all_25.gp
   PARI/GP's own verifier primecertisvalid applied to the 25 PARI/GP certificates the papers rely on
   (research/ergebnisse/w164_zertifikate/cert_1.txt ... cert_6.txt, Part I; w165_zertifikate/cert_1.txt ... cert_19.txt,
   Part II), so that these primes are checked by a second, independent design of checker besides ours.
   Run from the repository root:   gp -q checks/pari_all_25.gp
   Controls: PK: the control certificate cert_pk.txt (nextprime(10^40)) must be valid.
             NK: cert_1.txt of Part I with t of its first step shifted by 1 must be invalid.
   primecertisvalid returns 1 for a valid certificate and 0 otherwise; -1 below marks an error. */
default(parisizemax, 2000000000);
default(nbthreads, 1);
print("PARI/GP version: ", version());
check(label, f) =
{
  my(C = read(f), t0 = getwalltime(), r = iferr(primecertisvalid(C), E, -1));
  print(label, " ", f, ": ", #Str(C[1][1]), " digits, ", #C, " steps, ", getwalltime() - t0, " ms, result ", r);
  r == 1;
}
pk = check("PK", "research/ergebnisse/w164_zertifikate/cert_pk.txt");
nk = read("research/ergebnisse/w164_zertifikate/cert_1.txt");
nk[1][2] += 1;
nkr = iferr(primecertisvalid(nk), E, -1);
print("NK research/ergebnisse/w164_zertifikate/cert_1.txt with t of step 1 shifted by 1: result ", nkr);
ok = 0; n = 0;
for(i = 1, 6, n++; ok += check("Part I", Str("research/ergebnisse/w164_zertifikate/cert_", i, ".txt")));
for(i = 1, 19, n++; ok += check("Part II", Str("research/ergebnisse/w165_zertifikate/cert_", i, ".txt")));
print("RESULT: ", ok, " of ", n, " PARI/GP certificates accepted by primecertisvalid; positive control ", pk, ", negative control ", nkr);
