# SPDX-License-Identifier: Apache-2.0
"""Running times of the checker and the prover, written to a text file so that the numbers in README.md can be repeated.

    python tools/benchmark.py                 # certificates up to 700 digits, proofs up to 200 digits
    python tools/benchmark.py --all           # also the two PARI certificates of 1191 and 1687 digits (many minutes)

Times depend on the machine; the file records Python, platform, and whether gmpy2 is installed.
"""
import argparse
import hashlib
import pathlib
import platform
import random
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ecpp_research import arith, check_certificate, load  # noqa: E402
from ecpp_research.prove import prove  # noqa: E402

RESEARCH = ROOT / 'research'


def _key(path):
    """Folder first, then the number in the file name (cert_2 before cert_10; Part I before Part II)."""
    digits = ''.join(c if c.isdigit() else ' ' for c in path.stem).split()
    return (path.parent.name, 't2' in path.stem, int(digits[-1]) if digits else 0)


def certificates():
    pari = sorted((RESEARCH / 'ergebnisse').glob('w16[45]_zertifikate/cert_[0-9]*.txt'), key=_key)
    own = sorted((RESEARCH / 'ecpp' / 'ergebnisse').glob('e5_cert_*_eigen.txt'), key=_key)
    return [('PARI', p) for p in pari] + [('own', p) for p in own]


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--all', action='store_true', help='include certificates above 700 digits')
    parser.add_argument('--out', default=str(ROOT / 'tools' / 'benchmark_output.txt'))
    args = parser.parse_args(argv)
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    lines = []

    def say(text):
        print(text, flush=True)
        lines.append(text)

    say(f'Python {platform.python_version()} · {platform.platform()} · {platform.processor() or "cpu ?"} · '
        f'gmpy2 (prover only): {"yes" if arith.HAVE_GMPY2 else "no"}')
    say('')
    say('CHECKER (check_certificate): digits, steps, seconds')
    for origin, path in certificates():
        steps = load(path)
        digits = len(str(steps[0][0]))
        if digits > 700 and not args.all:
            say(f'  {origin:4s} {path.relative_to(ROOT).as_posix():58s} {digits:5d} digits  skipped (--all)')
            continue
        t0 = time.perf_counter()
        result = check_certificate(steps)
        say(f'  {origin:4s} {path.relative_to(ROOT).as_posix():58s} {digits:5d} digits {result.steps:4d} steps '
            f'{time.perf_counter() - t0:8.1f} s  {"VALID" if result.valid else "INVALID: " + result.reason}')
    say('')
    say('PROVER (prove, seed 1): N = next probable prime after 10^k; digits, steps, seconds, checked')
    rng = random.Random(1)
    for k in (40, 60, 100, 150, 200):
        n = 10 ** k + 1
        while not arith.probable_prime(n, rng):
            n += 2
        t0 = time.perf_counter()
        steps = prove(n, seed=1)
        seconds = time.perf_counter() - t0
        say(f'  10^{k} + {n - 10 ** k:<5d} {len(str(n)):5d} digits {len(steps):4d} steps {seconds:8.1f} s  '
            f'{"VALID" if check_certificate(steps, n).valid else "INVALID"}')
    text = '\n'.join(lines) + '\n'
    script = hashlib.md5(pathlib.Path(__file__).read_bytes()).hexdigest()[:8]
    text += f'md5 of this output (without this line) {hashlib.md5(text.encode("utf-8")).hexdigest()[:8]} · md5 of tools/benchmark.py {script}\n'
    pathlib.Path(args.out).write_text(text, encoding='utf-8', newline='\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
