\\ SPDX-License-Identifier: Apache-2.0
/* pari_minimal_reproduction.gp
   Six one-step certificates for composite N that PARI/GP's primecertisvalid accepts (2.17.2 and 2.19.0), reproduced
   with PARI/GP alone, without this package. For each certificate [N, t, s, a, [x, y]] it prints primecertisvalid, the
   factorization of N, q = (N + 1 - t)/s and isprime(q), and then, modulo each prime factor p of N, the order of the curve
   y^2 = x^3 + a*x + b (b from the point), s*P and q*(s*P). A valid step needs q*(s*P) to be the point at infinity
   (printed as [0]) modulo some prime factor; here it is not, modulo any of them. The last line is a positive control:
   a certificate made by primecert itself.
   Run:   gp -q checks/pari_minimal_reproduction.gp
   Found with checks/compare_with_pari.py (the case list) and reported to the PARI/GP developers on 25 October 2026 (bug #2723 in the PARI/GP bug tracker). */
print("version: ", version());
L = [[[5053, -74, 8, 4107, [965, 1506]]], [[4633, -67, 3, 3818, [1899, 3365]]], [[1921, 78, 4, 1678, [770, 1652]]], [[5587, 47, 3, 4416, [2461, 4582]]], [[10553, -194, 4, 1825, [4612, 1058]]], [[15931, 242, 30, 8239, [11912, 6118]]]];
{
for (i = 1, #L,
  C = L[i]; c = C[1]; N = c[1]; t = c[2]; s = c[3]; a = c[4]; x = c[5][1]; y = c[5][2];
  m = N + 1 - t; q = m / s; b = lift(Mod(y, N)^2 - Mod(x, N)^3 - a * Mod(x, N));
  print("certificate ", C, ": primecertisvalid = ", primecertisvalid(C), ", factor(N) = ", factor(N), ", q = ", q, ", isprime(q) = ", isprime(q));
  F = factor(N)[, 1];
  for (j = 1, #F,
    p = F[j]; E = ellinit([a, b], p); P = [x, y] * Mod(1, p);
    Q = ellmul(E, P, s); R = ellmul(E, Q, q);
    print("   mod ", p, ": #E = ", ellcard(E), ", s*P = ", lift(Q), ", q*(s*P) = ", lift(R))));
}
print("control, a certificate from primecert: ", primecertisvalid(primecert(nextprime(10^40))));
