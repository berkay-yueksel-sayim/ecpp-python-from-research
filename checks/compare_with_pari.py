# SPDX-License-Identifier: Apache-2.0
"""Compare the verdicts of this checker with PARI/GP's primecertisvalid on many small one-step certificates.

    python checks/compare_with_pari.py make                    # writes checks/compare_cases.gp and checks/compare_ours.txt
    gp -q checks/compare_with_pari.gp > checks/compare_pari.txt  # PARI/GP's verdicts (run from the repository root)
    python checks/compare_with_pari.py compare                 # compares, writes checks/compare_with_pari_output.txt

Cases (fixed seed): for prime N with 4 or 5 digits, the group order m = #E(F_N) is counted exactly, which gives valid
certificates [N, t, s, a, [x, y]] (t = N + 1 - m, q = m/s a prime above the bound, s*P not the point at infinity);
from each, invalid variants are made (t + 1, a + 1, y + 1, s replaced). For composite N = p*r, random certificates with
s dividing N + 1 - t are made, together with three certificates for composite numbers that the research checker and
the first draft of this package accepted. Expected: the same verdict from both checkers on every case, and no composite N
accepted by either.
"""
import math
import pathlib
import random
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ecpp_research import check_certificate, is_prime_small  # noqa: E402

HERE = ROOT / 'checks'
KNOWN_FALSE = [
    (10196873, 71, 1, 1021388, (4608684, 5613837)),
    (10210063, 3265, 1, 586389, (4233552, 698090)),
    (10270709, -1301, 1, 5725870, (4343939, 1328555)),
]


def legendre(v, p):
    v %= p
    return 0 if v == 0 else (1 if pow(v, (p - 1) // 2, p) == 1 else -1)


def group_order(a, b, p):
    return p + 1 + sum(legendre(x * x * x + a * x + b, p) for x in range(p))


def prime_factors(n):
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    if n > 1:
        out.append(n)
    return sorted(set(out))


def ref_mul(k, P, a, p):
    """k*P over F_p, exact group law."""
    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if (x1 - x2) % p == 0:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)
    R = None
    for bit in bin(k)[2:]:
        R = add(R, R)
        if bit == '1':
            R = add(R, P)
    return R


def bound(N):
    return (math.isqrt(math.isqrt(N)) + 2) ** 2


def make_cases(seed=20261023):
    rng = random.Random(seed)
    cases = []
    primes = [n for n in range(1000, 30000) if n % 6 in (1, 5) and is_prime_small(n)]
    while len(cases) < 1200:
        N = rng.choice(primes)
        a, x, y = rng.randrange(N), rng.randrange(N), rng.randrange(N)
        b = (y * y - x ** 3 - a * x) % N
        if (4 * a ** 3 + 27 * b * b) % N == 0:
            continue
        m = group_order(a, b, N)
        qs = [q for q in prime_factors(m) if q > bound(N)]
        if not qs:
            continue
        q = rng.choice(qs)
        s = m // q
        if ref_mul(s, (x, y), a, N) is None:
            continue
        t = N + 1 - m
        cases.append(('prime, valid', (N, t, s, a, (x, y))))
        for kind, step in (('prime, t+1', (N, t + 1, s, a, (x, y))), ('prime, a+1', (N, t, s, a + 1, (x, y))),
                           ('prime, y+1', (N, t, s, a, (x, y + 1))), ('prime, s*q', (N, t, s * q, a, (x, y)))):
            cases.append((kind, step))
    small = [p for p in range(5, 180) if is_prime_small(p)]
    while len(cases) < 2400:
        p, r = rng.sample(small, 2)
        N = p * r
        if N < 1000:
            continue
        t = rng.randrange(-math.isqrt(4 * N) + 1, math.isqrt(4 * N))
        m = N + 1 - t
        qs = [q for q in prime_factors(m) if q > bound(N)]
        if not qs:
            continue
        q = rng.choice(qs)
        cases.append(('composite, random', (N, t, m // q, rng.randrange(N), (rng.randrange(N), rng.randrange(N)))))
    for step in KNOWN_FALSE:
        cases.append(('composite, known false', step))
    return cases


def main(argv):
    if argv[:1] == ['make']:
        cases = make_cases()
        lines = [f'[[{N}, {t}, {s}, {a}, [{x}, {y}]]]' for kind, (N, t, s, a, (x, y)) in cases]
        (HERE / 'compare_cases.gp').write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')
        ours = [f'{kind}\t{int(check_certificate([step]).valid)}' for kind, step in cases]
        (HERE / 'compare_ours.txt').write_text('\n'.join(ours) + '\n', encoding='utf-8', newline='\n')
        print(f'{len(cases)} cases written')
        return 0
    if argv[:1] == ['compare']:
        ours = [z.split('\t') for z in (HERE / 'compare_ours.txt').read_text(encoding='utf-8').splitlines()]
        pari = [z.strip() for z in (HERE / 'compare_pari.txt').read_text(encoding='utf-8').splitlines()
                if z.strip() in ('-1', '0', '1')]  # -1: PARI/GP raised an error on the case (counts as a difference)
        assert len(pari) == len(ours), (len(pari), len(ours))
        table, differ = {}, []
        for i, ((kind, v), w) in enumerate(zip(ours, pari)):
            key = (kind, int(v), int(w))
            table[key] = table.get(key, 0) + 1
            if v != w:
                differ.append(i)
        ours_composite = sum(n for (kind, v, w), n in table.items() if kind.startswith('composite') and v == 1)
        pari_composite = sum(n for (kind, v, w), n in table.items() if kind.startswith('composite') and w == 1)
        lines = ['case kind | this checker | PARI/GP | count']
        lines += [f'{k} | {v} | {w} | {n}' for (k, v, w), n in sorted(table.items())]
        lines += ['', f'cases: {len(ours)}; verdicts differing: {len(differ)}; composite N accepted: this checker '
                      f'{ours_composite}, PARI/GP {pari_composite}', '']
        cases = make_cases()
        for i in differ:
            kind, (N, t, s, a, (x, y)) = cases[i]
            m, b = N + 1 - t, (y * y - x ** 3 - a * x) % N
            lines.append(f'case {i} ({kind}): [[{N}, {t}, {s}, {a}, [{x}, {y}]]], this checker {ours[i][1]}, PARI/GP {pari[i]}')
            if m % s == 0:
                q = m // s
                for p in prime_factors(N):
                    Q = ref_mul(s, (x % p, y % p), a % p, p)
                    R = ref_mul(q, Q, a % p, p) if Q is not None else None
                    lines.append(f'    modulo {p}: s*P = {Q}, q*(s*P) = {R}  (q = {q}; a valid step needs q*(s*P) = None, the point at infinity)')
        text = '\n'.join(lines) + '\n'
        (HERE / 'compare_with_pari_output.txt').write_text(text, encoding='utf-8', newline='\n')
        print(text)
        return 0 if ours_composite == 0 and all(not ours[i][0].startswith('prime') for i in differ) else 1
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
