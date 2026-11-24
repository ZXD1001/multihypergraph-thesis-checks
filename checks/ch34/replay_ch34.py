#!/usr/bin/env python3
"""Fixed finite checks for Chapter 3/4 examples; standard library only."""
from argparse import ArgumentParser
from collections import Counter
from decimal import Decimal, ROUND_CEILING, localcontext
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import platform
import sys


def inverse_mod(matrix, q):
    n = len(matrix)
    a = [list(row) + [int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j] % q)
        a[j], a[pivot] = a[pivot], a[j]
        inv = pow(a[j][j], -1, q)
        a[j] = [x * inv % q for x in a[j]]
        for i in range(n):
            if i != j:
                c = a[i][j]
                a[i] = [(x - c * y) % q for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


@lru_cache(maxsize=None)
def power_table(q):
    return [[pow(x, e, q) for e in range(q)] for x in range(q)]


def monomial(x, exponent, q):
    value = 1
    powers = power_table(q)
    for v, e in zip(x, exponent):
        value = value * powers[v][e] % q
    return value


def evaluate(x, coefficients, q):
    return sum(c * monomial(x, e, q) for e, c in coefficients.items()) % q


def recover_by_support(q, n, values):
    powers = power_table(q)
    b = [[powers[x][e] for e in range(1, q)] for x in range(1, q)]
    inverse = inverse_mod(b, q)
    recovered = {}
    block_log = []
    for t in range(1, n + 1):
        for support in combinations(range(n), t):
            support_set = set(support)
            tuples = list(product(range(1, q), repeat=t))
            residual = []
            for z in tuples:
                x = [0] * n
                for v, value in zip(support, z):
                    x[v] = value
                previous = {e: c for e, c in recovered.items()
                            if {i for i, a in enumerate(e) if a} < support_set}
                residual.append((values[tuple(x)] - evaluate(x, previous, q)) % q)
            block_coefficients = []
            for powers in tuples:
                value = 0
                for z, rhs in zip(tuples, residual):
                    entry = 1
                    for e, x in zip(powers, z):
                        entry = entry * inverse[e - 1][x - 1] % q
                    value += entry * rhs
                exponent = [0] * n
                for v, e in zip(support, powers):
                    exponent[v] = e
                recovered[tuple(exponent)] = value % q
                block_coefficients.append(value % q)
            block_log.append({'support': list(support), 'tuple_order': tuples,
                              'residual': residual, 'coefficients': block_coefficients})
    assert all(evaluate(x, recovered, q) == value for x, value in values.items())
    return recovered, block_log


def recovery_examples():
    cases = [
        ('q3_single_system', 1, {(1,): 2, (2,): 2}, 'eq:ch03-q3-one-variable'),
        ('q3_two_systems_delta', 2, {(1, 0): 2, (2, 0): 2, (1, 2): 1, (2, 2): 1},
         'app:w1-ordinary-hypergraph-constraints'),
        ('q3_three_systems_mg', 3, {(1, 1, 0): 1, (1, 2, 0): 2,
                                  (2, 2, 0): 2, (1, 0, 1): 1},
         'eq:ch03-fig32-mg-polynomial'),
        ('q3_three_systems_mhg', 3, {(1, 2, 0): 1, (1, 2, 1): 2,
                                   (2, 1, 2): 2, (0, 0, 2): 2},
         'eq:ch03-fig32-mhg-polynomial'),
    ]
    results = []
    for name, n, expected, label in cases:
        q = 3
        values = {x: evaluate(x, expected, q) for x in product(range(q), repeat=n)}
        assert values[(0,) * n] == 0
        if name == 'q3_two_systems_delta':
            assert all(value == int(x == (1, 0)) for x, value in values.items())
        recovered, blocks = recover_by_support(q, n, values)
        assert {e: c for e, c in recovered.items() if c} == expected
        results.append({'case': name, 'source_label': label, 'q': q, 'n': n,
                        'function_table': [{'x': x, 'value': v} for x, v in values.items()],
                        'recovered': [{'exponent': e, 'coefficient': c} for e, c in recovered.items()],
                        'blocks': blocks, 'full_back_substitution_points': len(values),
                        'expected_coefficients_match': True})
    return results


def ordinary_hg_effects():
    # Ordinary target F=x0*x1 on three qutrits. Compare all entries of
    # local measurement acceptance, the ZH relabeling rule and stabilizer averaging.
    q = 3
    points = list(product(range(q), repeat=3))
    reports = []
    for a in [(0, 2), (1,)]:
        b = tuple(v for v in range(3) if v not in a)
        checked = nonzero = 0
        for x in points:
            for y in points:
                checked += 1
                if any(x[v] != y[v] for v in b):
                    # All three expressions have zero entries here.
                    continue
                nonzero += 1
                cv = {0: x[1], 1: x[0], 2: 0}
                conditional_exponent = sum(cv[v] * (x[v] - y[v]) for v in a) % q
                accepted_fourier_labels = {v: -cv[v] % q for v in a}
                zh_exponent = sum(-accepted_fourier_labels[v] * (x[v] - y[v])
                                  for v in a) % q
                stabilizer_exponent = (x[0] * x[1] - y[0] * y[1]) % q
                assert conditional_exponent == zh_exponent == stabilizer_exponent
        reports.append({'A': a, 'B': b, 'checked_entries': checked,
                        'nonzero_entries': nonzero, 'common_denominator': q ** len(a),
                        'exact_root_exponents_equal_mod3': True})
    return {'q': 3, 'n': 3, 'target': 'F(x)=x0*x1', 'blocks': reports,
            'interpretation': 'fixed-target exact matrix-entry comparison; general equality uses the analytic reduction'}


def spectra_and_samples():
    rows = [('two_blocks_p_010', Fraction(1, 10), 2995, 2996),
            ('two_blocks_p_025', Fraction(1, 4), 1197, 1199),
            ('two_blocks_equal', Fraction(1, 2), 598, 600),
            ('three_singletons_equal', None, 898, 899)]
    results = []
    for name, weight, expected_n, expected_simple in rows:
        eigenvalues = {}
        for j in product(range(3), repeat=3):
            if weight is None:
                value = Fraction(sum(x == 0 for x in j), 3)
            else:
                value = weight * int(j[0] == j[2] == 0) + (1 - weight) * int(j[1] == 0)
            eigenvalues[j] = value
        gap = 1 - max(v for j, v in eigenvalues.items() if any(j))
        one_round_bound = 1 - gap * Fraction(1, 100)
        numerator, denominator = one_round_bound.numerator, one_round_bound.denominator
        assert 20 * pow(numerator, expected_n) <= pow(denominator, expected_n)
        assert 20 * pow(numerator, expected_n - 1) > pow(denominator, expected_n - 1)
        with localcontext() as ctx:
            ctx.prec = 80
            simple_real = Decimal(20).ln() * 100 * Decimal(gap.denominator) / Decimal(gap.numerator)
            simple = int(simple_real.to_integral_value(rounding=ROUND_CEILING))
        assert simple == expected_simple
        results.append({'case': name, 'p': str(weight) if weight is not None else None,
                        'spectrum_multiplicities': dict(Counter(str(v) for v in eigenvalues.values())),
                        'gap': str(gap), 'epsilon': '1/100', 'delta': '1/20',
                        'N': expected_n, 'N_and_Nminus1_checked_by_integer_powers': True,
                        'simplified_sufficient_N': simple,
                        'simplified_real_at_decimal_precision80': str(simple_real)})
    return results


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='New directory for detailed results and run metadata.')
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Do not use -O: assertions perform the verification.')
    args.output_dir.mkdir(parents=True, exist_ok=False)
    result = {'scope': 'four fixed parameter-recovery examples, two ordinary-HG blocks, four sample rows; not a universal implementation or proof',
              'recovery': recovery_examples(), 'ordinary_hg_effects': ordinary_hg_effects(),
              'spectra_samples': spectra_and_samples()}
    result_file = args.output_dir / 'ch34_checks.json'
    result_file.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    script = Path(__file__).resolve()
    manifest = {'python': sys.version, 'platform': platform.platform(),
                'script_sha256': sha256(script.read_bytes()).hexdigest(),
                'result_sha256': sha256(result_file.read_bytes()).hexdigest(),
                'passed': True, 'random_sampling': False}
    (args.output_dir/'run_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'recovery_cases': len(result['recovery']),
                      'back_substitution_points': sum(x['full_back_substitution_points'] for x in result['recovery']),
                      'ordinary_hg_matrix_entries': sum(x['checked_entries'] for x in result['ordinary_hg_effects']['blocks']),
                      'sample_rows': len(result['spectra_samples']), 'passed': True}, ensure_ascii=False))


if __name__ == '__main__':
    main()
