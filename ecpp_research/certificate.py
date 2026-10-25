# SPDX-License-Identifier: Apache-2.0
"""Reading and writing ECPP certificates in the PARI/GP ``primecert`` format.

A certificate is a list of steps ``[N, t, s, a, [x, y]]``. For step i the curve is
y^2 = x^3 + a*x + b with b = y^2 - x^3 - a*x (mod N), the group order is m = N + 1 - t, and
q = m / s is the number N of the next step. The last q is at most 2^64 and is proved prime by a
deterministic Miller-Rabin test. The format is the one described in the PARI/GP documentation
(``primecert``, ``primecertisvalid``); certificates written here can be checked by PARI/GP, and
PARI/GP certificates can be checked here.

Files hold the list as a literal: the same text is valid Python (``ast.literal_eval``) and GP.
"""
import ast
import contextlib
import pathlib
import sys


class CertificateFormatError(ValueError):
    """The text is not a list of steps [N, t, s, a, [x, y]] with integer entries."""


@contextlib.contextmanager
def whole_integers():
    """Lift Python's limit on converting long decimal integers (4300 digits by default since Python 3.11) for the
    duration of the block, and restore it afterwards. Certificates of numbers above 4300 digits need this."""
    get = getattr(sys, 'get_int_max_str_digits', None)
    if get is None:
        yield
        return
    old = get()
    sys.set_int_max_str_digits(0)
    try:
        yield
    finally:
        sys.set_int_max_str_digits(old)


def parse(text):
    """Parse certificate text and return a list of steps (N, t, s, a, (x, y)) with Python ints.

    Any text that is not such a list raises CertificateFormatError (also for malformed literals such as {[1]: 2},
    which ast.literal_eval rejects with TypeError, and for inputs too deeply nested or too large to evaluate)."""
    try:
        with whole_integers():
            data = ast.literal_eval(text.strip())
    except (ValueError, SyntaxError, TypeError, MemoryError, RecursionError) as exc:
        raise CertificateFormatError(f'not a certificate literal: {type(exc).__name__}: {exc}'[:300]) from None
    if not isinstance(data, (list, tuple)) or not data:
        raise CertificateFormatError('a certificate is a non-empty list of steps')
    steps = []
    for i, step in enumerate(data):
        if not isinstance(step, (list, tuple)) or len(step) != 5:
            raise CertificateFormatError(f'step {i}: expected [N, t, s, a, [x, y]]')
        N, t, s, a, P = step
        if not isinstance(P, (list, tuple)) or len(P) != 2:
            raise CertificateFormatError(f'step {i}: the point must be [x, y]')
        values = (N, t, s, a, P[0], P[1])
        if not all(isinstance(v, int) and not isinstance(v, bool) for v in values):
            raise CertificateFormatError(f'step {i}: all entries must be integers')
        steps.append((N, t, s, a, (P[0], P[1])))
    return steps


def load(path):
    """Read a certificate file."""
    return parse(pathlib.Path(path).read_text(encoding='utf-8'))


def dumps(steps):
    """Certificate text in the PARI/GP literal format (also valid Python)."""
    with whole_integers():
        return '[' + ', '.join(f'[{N}, {t}, {s}, {a}, [{P[0]}, {P[1]}]]' for N, t, s, a, P in steps) + ']'


def dump(steps, path):
    """Write a certificate file (UTF-8, newline at the end)."""
    pathlib.Path(path).write_text(dumps(steps) + '\n', encoding='utf-8', newline='\n')
