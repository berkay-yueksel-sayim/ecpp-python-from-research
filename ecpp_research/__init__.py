# SPDX-License-Identifier: Apache-2.0
"""ecpp_research: elliptic curve primality proving (ECPP) in pure Python, with an independent checker.

Written from the published descriptions of Goldwasser-Kilian (1986) and Atkin-Morain (1993) for a
paper series on consecutive powerful numbers; the research version and its certificates are kept
unchanged in ``research/``.

    from ecpp_research import check_certificate, load
    result = check_certificate(load('certificate.txt'))

    from ecpp_research.prove import prove          # needs mpmath
    steps = prove(2**127 - 1)                      # any probable prime above 2^64

The checker (``verify``, ``certificate``, ``curves``) needs only the standard library.
"""
from .certificate import CertificateFormatError, dump, dumps, load, parse
from .verify import Result, check_certificate, check_step, is_prime_small

__version__ = '0.1.1'
__all__ = ['CertificateFormatError', 'Result', 'check_certificate', 'check_step', 'dump', 'dumps',
           'is_prime_small', 'load', 'parse', '__version__']
