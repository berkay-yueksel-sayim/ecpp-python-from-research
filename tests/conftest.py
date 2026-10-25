# SPDX-License-Identifier: Apache-2.0
"""Shared test helpers: the repository root on sys.path, the research data, and the unchanged research checker w164."""
import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
RESEARCH = ROOT / 'research'
W164_FILE = RESEARCH / 'skripte' / 'w164_ecpp_pruefer_2026-10-02.py'


def load_research_checker():
    """Functions of the research checker w164, executed from its source text up to its main run (it has no import guard)."""
    src = W164_FILE.read_text(encoding='utf-8')
    head = src.split("pk = Z / 'cert_pk.txt'")[0]
    ns = {'__file__': str(W164_FILE), '__name__': 'w164_functions'}
    exec(compile(head, str(W164_FILE), 'exec'), ns)
    return ns


def research_certificates():
    """(label, path, expected first number or None) for every certificate kept in research/."""
    out = []
    z164 = RESEARCH / 'ergebnisse' / 'w164_zertifikate'
    out.append(('PARI nextprime(10^40)', z164 / 'cert_pk.txt', None))
    for job in json.loads((z164 / 'w164_auftrag.json').read_text(encoding='utf-8')):
        out.append((f'PARI Part I no. {job["i"]}', z164 / f'cert_{job["i"]}.txt', int(job['zahl'])))
    z165 = RESEARCH / 'ergebnisse' / 'w165_zertifikate'
    jobs2 = {j['i']: int(j['zahl']) for j in json.loads((z165 / 'w165_auftrag.json').read_text(encoding='utf-8'))}
    for i, n in sorted(jobs2.items()):
        out.append((f'PARI Part II no. {i}', z165 / f'cert_{i}.txt', n))
    own = RESEARCH / 'ecpp' / 'ergebnisse'
    jobs1 = {j['i']: int(j['zahl']) for j in json.loads((z164 / 'w164_auftrag.json').read_text(encoding='utf-8'))}
    out.append(('own nextprime(10^40)', own / 'e3_cert_pk_1e40.txt', None))
    for i, n in sorted(jobs1.items()):
        out.append((f'own Part I no. {i}', own / f'e5_cert_{i}_eigen.txt', n))
    for i, n in sorted(jobs2.items()):
        p = own / f'e5_cert_t2_{i}_eigen.txt'
        if p.exists():
            out.append((f'own Part II no. {i}', p, n))
    return out


@pytest.fixture(scope='session')
def w164():
    return load_research_checker()
