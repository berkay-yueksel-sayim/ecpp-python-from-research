# SPDX-License-Identifier: Apache-2.0
"""The checker: every certificate of the paper series is accepted, corruptions are rejected, and the result
agrees case by case with the unchanged research checker w164."""
import os
import random
import sys

import pytest

from conftest import research_certificates
from ecpp_research import check_certificate, check_step, dumps, is_prime_small, load, parse, CertificateFormatError

CERTS = research_certificates()
LARGE = 500  # certificates above this many digits run only with ECPP_SLOW=1 (each checked twice, up to minutes)


def _needs_slow(path):
    return len(str(load(path)[0][0])) > LARGE and not os.environ.get('ECPP_SLOW')


def test_inventory():
    """Expected inventory: 1 + 6 + 19 PARI certificates, 1 + 6 + 17 own ones (Part II nos. 6 and 16 have no own certificate)."""
    labels = [c[0] for c in CERTS]
    assert sum(l.startswith('PARI') for l in labels) == 26
    assert sum(l.startswith('own') for l in labels) == 24
    assert 'own Part II no. 6' not in labels and 'own Part II no. 16' not in labels


@pytest.mark.parametrize('label,path,number', CERTS, ids=[c[0] for c in CERTS])
def test_research_certificate_accepted(label, path, number, w164):
    if _needs_slow(path):
        pytest.skip('large certificate: set ECPP_SLOW=1')
    steps = load(path)
    result = check_certificate(steps, number)
    assert result.valid, result.reason
    assert w164['pruefe_zertifikat'](steps) is True  # the research checker agrees


def test_psi12_is_the_limit_of_the_small_test():
    """psi_12 = 399165290221 * 798330580441 (Sorenson-Webster, Thm 1.1) passes the 12-base test although composite;
    every q accepted at the end of a chain is <= 2^64 < psi_12."""
    psi12 = 318665857834031151167461
    assert psi12 == 399165290221 * 798330580441
    assert is_prime_small(psi12) is True
    assert psi12 > 2 ** 64
    assert is_prime_small(2 ** 61 - 1) and not is_prime_small(2 ** 61 + 1)


def test_corruptions_judged_like_w164(w164):
    """Random changes to the first step of small certificates.

    A change of N, t or s breaks the group order or the chain and must be rejected. A change of a or of the
    point need not: b is recomputed from the point, so the step then speaks about another curve, and for
    j = 0 or 1728 (a = 0 or b = 0, which PARI uses) a twist with the same group order is found with
    probability about 1/6 or 1/4 -- the changed certificate is then a valid proof again. For those fields
    only agreement with the research checker is required, and most changes must still be rejected.
    """
    rng = random.Random(20261021)
    small = [load(p) for label, p, n in CERTS if len(str(load(p)[0][0])) <= 120]
    assert len(small) >= 5
    rejected = total = 0
    for steps in small:
        for _ in range(40):
            bad = [list(s) for s in steps]
            field = rng.randrange(6)
            if field == 5:
                bad[0][4] = (bad[0][4][0], bad[0][4][1] + rng.randrange(1, 10))
            elif field == 4:
                bad[0][4] = (bad[0][4][0] + rng.randrange(1, 10), bad[0][4][1])
            elif field == 0:
                bad[0][0] += 2
            else:
                bad[0][field] += rng.choice([-2, -1, 1, 2])
            bad = [tuple(s) for s in bad]
            new = check_certificate(bad).valid
            assert new == (w164['pruefe_zertifikat'](bad) is True)
            if field in (0, 1, 2):
                assert not new
            total += 1
            rejected += not new
    assert rejected >= 0.8 * total


def test_step_agrees_with_w164_on_every_step(w164):
    """For every step of the certificates up to 250 digits, both checkers return the same q."""
    for label, p, n in CERTS:
        steps = load(p)
        if len(str(steps[0][0])) > 250:
            continue
        for N, t, s, a, P in steps:
            assert check_step(N, t, s, a, P) == w164['pruefe_schritt'](N, t, s, a, P)


def test_chain_must_link():
    steps = load(CERTS[0][1])
    broken = [steps[0], steps[2]] + steps[3:] if len(steps) > 2 else steps
    assert not check_certificate(broken).valid


def test_format_roundtrip_and_errors():
    steps = load(CERTS[0][1])
    assert parse(dumps(steps)) == steps
    for bad in ['', '[]', '[[1, 2, 3]]', '[[5, 1, 1, 1, [True, 2]]]', 'import os', '{[1]: 2}', '[' * 100000]:
        with pytest.raises(CertificateFormatError):
            parse(bad)


def test_more_than_4300_digits():
    """Python refuses decimal integers above 4300 digits by default (3.11+); reading and writing lift the limit for
    the duration of the call only."""
    big = 10 ** 5000 + 1
    steps = [(big, 0, 7, 1, (2, 3))]  # m = 10^5000 + 2 = 4 (mod 7), so s = 7 does not divide m
    assert parse(dumps(steps)) == steps
    assert check_certificate(steps).reason.endswith('s does not divide m')  # rejected, not a crash
    if hasattr(sys, 'get_int_max_str_digits'):
        before = sys.get_int_max_str_digits()
        parse(dumps(steps))
        assert sys.get_int_max_str_digits() == before
