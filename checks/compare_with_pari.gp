\\ SPDX-License-Identifier: Apache-2.0
/* compare_with_pari.gp
   PARI/GP's verdict (1 valid, 0 invalid) for every certificate in checks/compare_cases.gp, one per line, in order.
   Run from the repository root:   gp -q checks/compare_with_pari.gp > checks/compare_pari.txt
   The cases and this checker's verdicts are made by:   python checks/compare_with_pari.py make */
L = readvec("checks/compare_cases.gp");
for(i = 1, #L, print(iferr(primecertisvalid(L[i]), E, -1)));
