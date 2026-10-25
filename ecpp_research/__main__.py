# SPDX-License-Identifier: Apache-2.0
"""Command line:

    python -m ecpp_research verify CERTIFICATE [--number N]
    python -m ecpp_research prove N [--out FILE] [--time-limit SECONDS] [--seed S] [-v]

``verify`` needs only the Python standard library. ``prove`` needs mpmath (and uses gmpy2 if installed).
Exit codes: 0 the certificate is valid / the proof was found and verified; 1 the certificate is invalid / no proof
(N <= 2^64, not a probable prime, or no certificate within the limits); 2 the input cannot be read or is not an integer.
"""
import argparse
import sys
import time

from .certificate import CertificateFormatError, dump, load, whole_integers
from .verify import check_certificate


def _verify(args):
    try:
        steps = load(args.certificate)
    except (OSError, CertificateFormatError) as exc:
        print(f'cannot read certificate: {exc}', file=sys.stderr)
        return 2
    t0 = time.time()
    result = check_certificate(steps, args.number)
    text = str(result.number)
    shown = text if len(text) <= 40 else f'{text[:20]}...{text[-12:]}'
    if result.valid:
        print(f'VALID: {shown} ({len(text)} digits) proved prime ({result.steps} steps, {time.time() - t0:.1f} s)')
        return 0
    print(f'INVALID: {result.reason}')
    return 1


def _prove(args):
    from .prove import ProofNotFound, prove  # imports mpmath only when proving
    try:
        n = int(args.number)
    except ValueError:
        print(f'not an integer: {args.number[:60]}', file=sys.stderr)
        return 2
    log = (lambda line: print('  ' + line, file=sys.stderr, flush=True)) if args.verbose else None
    t0 = time.time()
    try:
        steps = prove(n, time_limit=args.time_limit, seed=args.seed, log=log)
    except ValueError as exc:
        print(f'not proved: {exc}')
        return 1
    except ProofNotFound as exc:
        print(f'no proof found: {exc}')
        return 1
    if args.out:
        dump(steps, args.out)
    print(f'PROVED: {len(str(n))}-digit number, {len(steps)} steps, {time.time() - t0:.1f} s'
          + (f', certificate written to {args.out}' if args.out else ''))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog='ecpp-research', description='ECPP prover and independent certificate checker.')
    sub = parser.add_subparsers(dest='command', required=True)
    v = sub.add_parser('verify', help='check a certificate file (PARI/GP primecert format)')
    v.add_argument('certificate')
    v.add_argument('--number', type=int, default=None, help='require that the certificate proves this number')
    v.set_defaults(func=_verify)
    p = sub.add_parser('prove', help='build a certificate for a number > 2^64')
    p.add_argument('number')
    p.add_argument('--out', default=None, help='write the certificate to this file')
    p.add_argument('--time-limit', type=float, default=3600.0)
    p.add_argument('--seed', type=int, default=None)
    p.add_argument('-v', '--verbose', action='store_true', help='print progress to stderr')
    p.set_defaults(func=_prove)
    with whole_integers():  # numbers above 4300 digits on the command line and in the printed digit counts
        args = parser.parse_args(argv)
        return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
