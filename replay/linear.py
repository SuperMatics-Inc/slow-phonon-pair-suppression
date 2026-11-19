"""Exact rational subchecks; --proxies operates on the inherited decimal table."""

import argparse
import csv
from fractions import Fraction as Q
from pathlib import Path

MUTANTS = (
    'frequency-sign', 'endpoint-count', 'attachment', 'orientation',
    'temperature-power', 'root-split', 'strip-cutoff', 'dos-per-spin',
)
LH, LS, KH, KS = 2**26, 2**40, 2**34, 2**74


def check(condition, reason):
    if not condition:
        raise SystemExit(f'REJECT: {reason}')


def distance(left, right):
    return sum(abs(x - y) for x, y in zip(left, right))


def replay(mutant):
    odd_signs = 1 if mutant == 'frequency-sign' else 2
    basel_over_pi_squared = Q(1, 6)
    inverse_frequency_mass = odd_signs * (1 - Q(1, 4)) * basel_over_pi_squared
    check(inverse_frequency_mass == Q(1, 4), 'two signs in the odd inverse-square sum')
    weighted_operator = 1 + 2 * Q(1, 2) * inverse_frequency_mass
    check(weighted_operator == Q(5, 4), 'weighted operator norm')

    shifts = 1 if mutant == 'endpoint-count' else 2
    left, right, common_leg = (1, 0), (0, 1), (1, 1)
    check(distance(left, right) <= distance(common_leg, common_leg) + shifts,
          'both anchor offsets can be nonzero')
    envelope = 16 if mutant == 'attachment' else 64
    check(3 * weighted_operator * 8 <= envelope and 3 <= 8 and 32 == 8 * 4,
          'weighted row sum from the two-dimensional Cauchy–Schwarz estimate')
    check(1 + 4 * 4097 <= 2**15, 'reference hopping budget')
    ssh_orientations = 1 if mutant == 'orientation' else 2
    check(envelope * 4 * 4096 * 64 == LH
          and envelope * 2**15 * 4096 * ssh_orientations * 64 == LS,
          'both quantum budgets, including the SSH orientation factor')
    check(envelope * 2**15 * 8 == 2**24 and 2**21 - 1 > 24 + 6,
          'small uniform SSH tail')

    power = 3 if mutant == 'temperature-power' else 4
    for tau in (Q(1, 3), Q(1), Q(11, 5)):
        for energy in (Q(1, 2), Q(17, 4)):
            beta = 1 / tau
            check(beta**3 * (1 + beta * energy) * tau**power == energy + tau,
                  'conversion of four inverse powers of temperature')
    budget = Q(2, 3) if mutant == 'root-split' else Q(1, 2)
    check(2 * budget <= 1, 'two terms of the explicit root upper bound')

    ylo, yhi, elo, ehi = Q(1, 4), Q(1, 2), 1, 2
    cutoff = 32 if mutant == 'strip-cutoff' else 64
    check(ylo**2 / (ehi + yhi) == Q(1, 40) > Q(1, cutoff), 'new strip lies above sharpened test')
    for old, new in ((KH, LH), (KS, LS)):
        check(yhi**4 / elo**3 == Q(1, 16) < Q(old, 2 * cutoff * new),
              'new strip lies below old test')
        check(Q(old, 128 * new) >= 2 and 3 * Q(3, 4) / 4 == Q(9, 16) > yhi,
              'old bound and G3 miss the strip')
        check(ylo**2 / (128 * new) == Q(1, 2048 * new)
              and yhi**2 / (64 * new * Q(1, 2)) == Q(1, 128 * new),
              'frequency window of the strip')

    # Jacobian of E=4tw, spin factor for n, and power-majorant integral at w0.
    spin = 1 if mutant == 'dos-per-spin' else 2
    check(4 * spin == 8 and 8 * Q(1, 4) * 2 == 4, 'per-spin DOS to total density')
    sqrt_w = Q(1, 32)
    check(sqrt_w**2 == Q(1, 1024) and 4 * sqrt_w / 3 < Q(1, 4)
          and 4 * 2 / (4 * 3 * sqrt_w) < 24 < 64 * LH <= 64 * LS,
          'a regular-density endpoint below the open strip')
    print('PASS: weighted covariance, temperature, strip and DOS arithmetic')


def proxies(csv_path):
    with csv_path.open(newline='', encoding='utf-8') as stream:
        rows = [row for row in csv.DictReader(stream) if row['source'] == 'PRX Fig. 9']
    cinf = Q('0.18272624777232952')  # used only to obtain the simple bound C_inf<1/5
    check(cinf < Q(1, 5), 'bare comparator is less than 1/5')
    print('NUMERICAL (inherited digitized inputs): omega/t, Hubbard/t, lambda_FS, a_S')
    for omega in (Q(3, 10), Q(1, 2), Q(1)):
        for hubbard in (Q(0), Q(8)):
            candidates = [row for row in rows if Q(row['omega']) == omega and Q(row['u']) == hubbard]
            row = max(candidates, key=lambda point: Q(point['r_wide_184']))
            lam_fs = Q(row['moment_b']) / omega**2
            a_s = lam_fs / (4 * Q(row['lam']))
            # Y^2 < omega^2*lambda_FS/25 <= 2/25 and its divisor is >= 1.
            # The RHS of (2) exceeds LS when a_S < 1.
            check(0 < lam_fs < 2 and 0 < a_s < 1, 'six candidate rows are outside (2)')
            print(f'{float(omega):g}, {float(hubbard):g}, {float(lam_fs):.6g}, {float(a_s):.6g}')


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--mutant', choices=MUTANTS)
    cli.add_argument('--proxies', type=Path, metavar='PUBLISHED_MAP_CSV')
    arguments = cli.parse_args()
    replay(arguments.mutant)
    if arguments.proxies:
        proxies(arguments.proxies)
