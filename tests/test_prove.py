# SPDX-License-Identifier: Apache-2.0
"""Controls of the prover, carried over from the research script e3 (K1-K5), plus the cross route:
every certificate the package builds must also be accepted by the unchanged research checker w164."""
import os
import random

import pytest

from conftest import ROOT
from ecpp_research import arith, check_certificate, check_step, load
from ecpp_research.__main__ import main
from ecpp_research.prove import ProofNotFound, prove, step
from ecpp_research.classpoly import discriminants


def next_probable_prime(n, rng):
    n += 1 - n % 2
    while not arith.probable_prime(n, rng):
        n += 2
    return n


def accepted_by_both(steps, number, w164):
    ours = check_certificate(steps, number)
    research = w164['pruefe_zertifikat']([tuple(s[:4]) + (list(s[4]),) for s in steps])
    return ours.valid and research is True


def test_reference_point_counts():  # e3 K5 (Atkin-Morain pp. 59-60)
    for p, D, count in ((17401, 8, 17200), (107, 7, 88)):
        A, B = arith.cornacchia4(p, D)
        assert count in (p + 1 - A, p + 1 + A)


def test_nextprime_1e40(w164):  # e3 K1
    rng = random.Random(40)
    N = next_probable_prime(10 ** 40, rng)
    pari = load(ROOT / 'research' / 'ergebnisse' / 'w164_zertifikate' / 'cert_pk.txt')
    assert pari[0][0] == N  # the number PARI's control certificate proves
    steps = prove(N, time_limit=300, seed=40)
    assert accepted_by_both(steps, N, w164)


@pytest.mark.parametrize('digits', [30, 45, 60])
def test_random_primes(digits, w164):  # e3 K2
    rng = random.Random(digits)
    for _ in range(3):
        N = next_probable_prime(rng.randrange(10 ** (digits - 1), 10 ** digits), rng)
        steps = prove(N, time_limit=300, seed=rng.randrange(2 ** 32))
        assert accepted_by_both(steps, N, w164)


def test_composites():  # e3 K3
    rng = random.Random(3)
    D = discriminants(8, 6000)[:12]
    for _ in range(3):
        n = (next_probable_prime(rng.randrange(10 ** 19, 10 ** 20), rng)
             * next_probable_prime(rng.randrange(10 ** 19, 10 ** 20), rng))
        assert step(n, rng, D) is None
        with pytest.raises(ValueError):
            prove(n)
    with pytest.raises(ValueError):
        prove(2 ** 61 - 1)  # below 2^64: the deterministic small test is the right tool there


def test_falsified_first_step(w164):  # e3 K4
    rng = random.Random(4)
    N = next_probable_prime(10 ** 40, rng)
    first = list(prove(N, time_limit=300, seed=4)[0])
    variants = [list(first) for _ in range(3)]
    variants[0][1] += 1                                  # t shifted by 1
    variants[1][2] *= 2                                  # s doubled
    variants[2][4] = (first[4][0], first[4][1] + 1)      # point moved off the curve
    for v in variants:
        assert isinstance(check_step(*v[:4], v[4]), str)
        assert isinstance(w164['pruefe_schritt'](*v[:4], list(v[4])), str)


def test_same_seed_same_certificate():
    rng = random.Random(5)
    N = next_probable_prime(10 ** 35, rng)
    assert prove(N, seed=11) == prove(N, seed=11)


def test_time_limit_is_reported():
    rng = random.Random(6)
    N = next_probable_prime(10 ** 60, rng)
    with pytest.raises(ProofNotFound):
        prove(N, time_limit=0.0)


def test_command_line(tmp_path, capsys):
    rng = random.Random(7)
    N = next_probable_prime(10 ** 30, rng)
    cert = tmp_path / 'cert.txt'
    assert main(['prove', str(N), '--out', str(cert), '--seed', '1']) == 0
    assert main(['verify', str(cert), '--number', str(N)]) == 0
    assert main(['verify', str(cert), '--number', str(N + 2)]) == 1
    cert.write_text('not a certificate', encoding='utf-8')
    assert main(['verify', str(cert)]) == 2
    assert main(['prove', str(N + 1)]) == 1  # even, hence not a probable prime
    assert 'VALID' in capsys.readouterr().out


@pytest.mark.skipif(not os.environ.get('ECPP_SLOW'), reason='minutes: set ECPP_SLOW=1')
def test_part_one_numbers_reproved(w164):
    """Integration: re-prove the six numbers of Part I (111-477 digits) from scratch."""
    folder = ROOT / 'research' / 'ecpp' / 'ergebnisse'
    for k in range(1, 7):
        N = load(folder / f'e5_cert_{k}_eigen.txt')[0][0]
        steps = prove(N, time_limit=3600, seed=k)
        assert accepted_by_both(steps, N, w164)
