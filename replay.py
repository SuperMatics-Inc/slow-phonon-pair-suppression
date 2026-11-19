"""Exact finite scalar checks for the consolidated manuscript (stdlib only)."""

import argparse
from fractions import Fraction as Q
from pathlib import Path
import runpy
import subprocess
import sys

HERE = Path(__file__).resolve().parent
MUTANTS = ('spectral-diagonal', 'spectral-upper', 'ivt-density', 'dyadic-exponent',
           'dos-factor', 'ssh-stiffness')


def require(condition, guard):
    if not condition:
        raise SystemExit(f'REJECT: {guard}')


def spectral_bracket(mutant):
    pi_lo, pi_hi = Q(157, 50), Q(22, 7)
    c_lo, c_hi = Q(9, 50), Q(1, 5)
    v = (Q(1), Q(2, 3), Q(1, 4))
    sign = 1 if mutant == 'spectral-diagonal' else -1

    def numerator(i, j):
        entry = Q(1, (i + j + 1)**2)
        if i != j:
            return entry + Q(1, (i - j)**2)
        return entry + sign * 2 * sum(Q(1, k*k) for k in range(1, i + 1))

    require(numerator(1, 1) / 3 == -Q(17, 27), 'soft-kernel normal diagonal subtraction')
    rayleigh = (sum(v[i]*v[j]*numerator(i, j)/((2*i+1)*(2*j+1))
                    for i in range(3) for j in range(3))
                 / sum(v[i]**2/(2*i+1) for i in range(3)))
    require(rayleigh > (2*pi_hi*c_lo)**2, 'spectral trial lower bound')
    require(sum(Q(1, k*k) for k in range(1, 6)) + Q(1, 5) < Q(5, 3), 'inverse-square sum bound')
    eta0 = Q(1) if mutant == 'spectral-upper' else Q(5, 4)
    rows = (1 + (Q(5, 12) + (Q(3, 4) + Q(4, 9))/5) / eta0,
            -Q(17, 27) + eta0*Q(5, 4) + Q(2, 5),
            eta0*(Q(1, 4) + Q(1, 9)) + Q(17, 48) + Q(2, 3) + Q(1, 20))
    require(max(rows) < (2*pi_lo*c_hi)**2, 'signed Schur upper bound')
    return c_lo, c_hi


def witnesses(mutant):
    c_lo, c_hi = spectral_bracket(mutant)
    u = Q(3, 4)
    ivt_target = 48 if mutant == 'ivt-density' else 96
    exponent = 36000 if mutant == 'dyadic-exponent' else 72000
    spin = 1 if mutant == 'dos-factor' else 2
    require(4*spin == 8, 'per-spin DOS contributes both spins to n')
    wmax = Q(1, 1024)
    require(4*Q(1, 32)/3 < Q(1, 4) and 2*(1+wmax) + 4*wmax < 4,
            'filling endpoint and upper logarithmic correction')
    require(Q(2, 3)/(3**3) / (1-Q(1, 9)) + Q(2, 3) == Q(25, 36) < Q(7, 10),
            'log 2 geometric-series enclosure')
    require(Q(99, 1000) < 1/(Q(22, 7)**2*(1+wmax))
            and (1+wmax*wmax)/9 < Q(101, 900), 'DOS prefactor envelope')
    require((c_lo**2 * Q(99, 1000) * (24000-8), c_hi**2 * Q(101, 900) * (25200+8))
            == (Q(24048981, 312500), Q(636502, 5625)), 'exact dyadic product endpoints')
    for lconst, bond_factor in ((2**26, 1), (2**40, 4)):
        require(c_hi**2 * 24 < ivt_target*lconst and 64 < ivt_target < 128
                and ivt_target*lconst*u / (512*lconst) == Q(3, 8)**2, 'IVT witness')
        require(exponent*lconst % bond_factor == 0 and exponent*lconst//bond_factor > 10,
                'dyadic exponent is an integer away from endpoint')
        lo = c_lo**2 * Q(99, 1000) * (Q(exponent, 3)*lconst - 8)
        hi = c_hi**2 * Q(101, 900) * (Q(7*exponent, 20)*lconst + 8)
        require(64*lconst < lo < hi < 128*lconst, 'dyadic witness in strict strip')
    ssh_multiplier = 1+u if mutant == 'ssh-stiffness' else 1+2*u
    require(Q(1, 2)*ssh_multiplier == Q(1, 2)+u, 'SSH phase stiffness has the full phonon term')
    require(Q(1, 4)**2/(2+Q(1, 2)) > Q(1, 64) and 3*Q(3, 4)/4 > Q(1, 2),
            'strict adiabatic comparison and both stiffness crossings')
    print('PASS: rational spectral bracket, IVT and dyadic witnesses, SSH multiplier')


def run_suite(run_mutants):
    for name, path, mutants in (
        ('base', HERE/'replay'/'base.py', None),
        ('linear', HERE/'replay'/'linear.py', None),
        ('witness', HERE/'replay.py', MUTANTS),
    ):
        if mutants is None:
            mutants = runpy.run_path(str(path))['MUTANTS']
            subprocess.run([sys.executable, '-B', str(path)], check=True)
        if run_mutants:
            for mutant in mutants:
                result = subprocess.run([sys.executable, '-B', str(path), '--mutant', mutant],
                                        capture_output=True, text=True)
                require(result.returncode == 1 and result.stderr.startswith('REJECT:'),
                        f'{name}/{mutant} did not hit an intended guard: {result.stderr}')
                print(f'{name}/{mutant} exit 1: {result.stderr.strip()}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--all-mutants', action='store_true')
    modes.add_argument('--mutant', choices=MUTANTS)
    args = parser.parse_args()
    witnesses(args.mutant)
    if args.mutant is None:
        run_suite(args.all_mutants)
