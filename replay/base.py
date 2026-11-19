"""Exact replay of finite ADIABATIC guards; proxy table is explicitly floating point."""
import argparse
import csv
from fractions import Fraction as F
from itertools import combinations, permutations
from math import pi, sqrt
from pathlib import Path


MUTANTS = ('phonon-sign', 'matsubara-pair', 'gaussian-half', 'ssh-orientation',
           'cutoff-tail', 'attachment', 'endpoint-count', 'pair-charge', 'jump-factor',
           'quadratic', 'forest-path')


def require(condition, guard):
    if not condition:
        raise SystemExit(f'REJECT: {guard}')


def forest_monomial(edges, mutant):
    """All listed edges are distinct copies of a connected multigraph."""
    vertices = {x for edge in edges for x in edge}
    total = F(0)
    for tree in combinations(range(len(edges)), len(vertices) - 1):
        adj = {x: [] for x in vertices}
        for i in tree:
            x, y = edges[i]
            adj[x].append((y, i))
            adj[y].append((x, i))

        def path(start, end, seen=frozenset()):
            if start == end:
                return ()
            for y, edge in adj[start]:
                if y not in seen:
                    suffix = path(y, end, seen | {start})
                    if suffix is not None:
                        return (edge,) + suffix
            return None

        if any(path(min(vertices), x) is None for x in vertices):
            continue
        paths = [path(*edge) for i, edge in enumerate(edges) if i not in tree]
        for order in permutations(tree):
            powers = {i: 0 for i in tree}
            for p in paths:
                powers[next(i for i in (order[::-1] if mutant else order) if i in p)] += 1
            accumulated, volume = 0, F(1)
            for rank, i in enumerate(order, 1):
                accumulated += powers[i]
                volume /= rank + accumulated
            total += volume
    return total


def verify(mutant):
    # Coefficient of U β Ω² after the exact Basel sum and ±k doubling.
    sign = -1 if mutant == 'phonon-sign' else 1
    pairs = 1 if mutant == 'matsubara-pair' else 2
    dstar = sign * pairs * F(1, 4) * F(1, 6)
    require(dstar == F(1, 12), 'Fourier-pair covariance coefficient')
    half = F(1) if mutant == 'gaussian-half' else F(1, 2)
    require(half * 4**2 * dstar == F(2, 3), 'exact Gaussian vertex coefficient')
    e_upper = F(11, 4)  # certified below by a geometric tail
    partial, factorial = F(1), 1
    for k in range(1, 6):
        factorial *= k
        partial += F(1, factorial)
    require(2 < partial and partial + F(1, 720) * F(7, 6) < e_upper, 'Euler-number bracket')
    require(4 * e_upper**4 < 256 and F(2, 3)*e_upper**8 < 4096, 'monomial norm bounds')
    orientations = 1 if mutant == 'ssh-orientation' else 2
    require(orientations == 2, 'both bond orientations')
    decay_sign = -1 if mutant == 'cutoff-tail' else 1
    require(decay_sign * 4096**2 // 8 == 2**21 and 4096 >= 4*512,
            'Gaussian-tail exponent and absorbing radius')
    endpoint_shifts = 1 if mutant == 'endpoint-count' else 2
    left_anchor, right_anchor, shared_end = (1, 0), (0, 1), (1, 1)
    def distance(x, y):
        return sum(abs(a-b) for a, b in zip(x, y))
    require(distance(left_anchor, right_anchor) <= distance(shared_end, shared_end)
            + endpoint_shifts, 'both SSH endpoint shifts')
    envelope = 512 if mutant == 'attachment' else 1024
    require(e_upper < 3 and e_upper**2 < 8 and F(36, 5) < 8
            and max(3 * 16**2, 8 * 8**2) <= envelope,
            'two-case weighted square lattice attachment budget')
    require(4*4 <= 16 and 4 <= 16*3 and 36-5*16 + F(1, 2)*(36*3-5*4) == 0
            and 36*3-5*4 >= 0, 'rational attachment ratio in both intervals')
    require(1 + 4*4097 <= 2**15, 'SSH cutoff hopping budget')
    q_h, q_s = envelope * 4**3, envelope * (2**15)**3
    require(q_h == 2**16 and q_s == 2**55, 'attachment constants')
    require(q_h * 4096 * 64 == 2**34, 'Holstein temperature constant')
    require(q_s * 4096 * orientations * 64 == 2**74, 'SSH temperature constant')
    require(q_s * 8 == 2**58 and 2**21 - 3 > 58 + 6, 'SSH strict Gaussian-tail margin')
    require(F(1, 64) + F(1, 64) < F(1, 4), 'additional small-source margin')
    charge = 1 if mutant == 'pair-charge' else 2
    jump = F(1) if mutant == 'jump-factor' else F(1, 2)  # in units of π
    require(F(1, charge**2) == F(1, 4), 'pair gradient / electron Peierls factor')
    require(jump * F(2, charge**2) == F(1, 4), 'G3 BKT comparison prefactor')
    c_upper, a_s, y_limit = F(1, 5), F(14, 25), F(13)
    require(c_upper * 32 * 2 < y_limit < F(2**34) / (c_upper**2 * a_s),
            'dynamic linear comparison off-scale in both illustrated panels')
    for q, e in ((F(1), F(0)), (F(1), F(2)), (F(1, 9), F(10, 9))):
        # Exact polynomial coefficients of (2τ-q)^2 - (q^2+4qE).
        b = 4 if mutant == 'quadratic' else -4
        require((4, b*q, -4*q*e) == tuple(4*x for x in (1, -q, -q*e)),
                'positive quadratic root discriminant / sign')
    for graph in (((0, 1), (0, 1)), ((0, 1), (1, 2), (0, 2)),
                  ((0, 1), (1, 2), (0, 2), (0, 1), (1, 3), (2, 3))):
        require(forest_monomial(graph, mutant == 'forest-path') == 1,
                'maximum-spanning-tree polynomial identity')
    print('PASS: exact term, vertex, attachment, tail, temperature, flux and forest guards')


def proxies(path):
    with path.open(newline='', encoding='utf-8') as file:
        rows = list(csv.DictReader(file))
    cinf = 0.18272624777232952
    print('NUMERICAL, inherited digitized inputs: omega/t, lambda_Z, Hubbard/t, n, '
          'lambda_FS, G3 crossing omega/t, central proxy ratio')
    for omega in (0.3, 0.5, 1.0):
        for hubbard in (0.0, 8.0):
            choices = [r for r in rows if r['source'] == 'PRX Fig. 9'
                       and float(r['omega']) == omega and float(r['u']) == hubbard]
            row = max(choices, key=lambda r: float(r['r_wide_184']))
            lam_z, n, lambda_fs = (float(row[k]) for k in ('lam', 'density_b', 'moment_b'))
            lambda_fs /= omega**2
            oscillator_u = 4 * lam_z
            crossing = pi * (1 + 2 * oscillator_u) * min(n, 2-n) / (4*cinf*sqrt(lambda_fs))
            print(*(f'{x:.5g}' for x in (omega, lam_z, hubbard, n, lambda_fs,
                                         crossing, float(row['r_b']))), sep=', ')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mutant', choices=MUTANTS)
    parser.add_argument('--proxies', type=Path, metavar='PUBLISHED_MAP_CSV')
    args = parser.parse_args()
    verify(args.mutant)
    if args.proxies:
        proxies(args.proxies)
