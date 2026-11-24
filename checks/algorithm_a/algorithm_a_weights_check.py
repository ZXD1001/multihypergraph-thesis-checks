"""Bounded exact checks for the two interpolation weights and D_i S_i order.

Two finite examples check Algorithm A interpolation and the order of operations.
The enumeration does not establish correctness for arbitrary parameters.
Only the Python standard library is required.
"""
from itertools import product
from math import prod
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys


def check_case(q, k, m, points, shift, high_functions, linear, quadratic):
    n = 2 * k - 1
    mask_count = m * (k - 1)
    assert len(points) == n and len(set(points)) == n
    assert len(shift) == mask_count == len(high_functions)
    nu = [
        pow(prod((points[i] - points[j]) % q for j in range(k) if j != i) % q, -1, q)
        if i < k else 0
        for i in range(n)
    ]
    lam = [
        pow(prod((points[i] - points[j]) % q for j in range(n) if j != i) % q, -1, q)
        for i in range(n)
    ]
    assert nu != lam  # The example must distinguish the two sets of weights.
    assert all(sum(nu[i] * pow(points[i], d, q) for i in range(n)) % q == int(d == k - 1)
               for d in range(k))
    assert all(sum(lam[i] * pow(points[i], d, q) for i in range(n)) % q == int(d == 2*k - 2)
               for d in range(2*k - 1))

    def values(mask, secret):
        return [[(sum(mask[r*(k-1)+j] * pow(alpha, j, q) for j in range(k-1))
                  + secret[r] * pow(alpha, k-1, q)) % q
                 for r in range(m)] for alpha in points]

    def h(mask, secret):
        return sum(a * f(secret) for a, f in zip(mask, high_functions)) % q

    def low(secret):
        return (sum(g*x for g, x in zip(linear, secret))
                + sum(g*secret[r]*secret[s] for r, s, g in quadratic)) % q

    delta = [[sum(shift[r*(k-1)+j]*pow(alpha, j, q) for j in range(k-1)) % q
              for r in range(m)] for alpha in points]
    count = 0
    for coords in product(range(q), repeat=mask_count+m):
        mask, secret = coords[:mask_count], coords[mask_count:]
        before_mask = tuple((a+t) % q for a, t in zip(mask, shift))
        post_values = values(mask, secret)
        pre_values = values(before_mask, secret)
        assert all((pre_values[i][r] - delta[i][r]) % q == post_values[i][r]
                   for i in range(n) for r in range(m))
        physical_phase = sum(
            nu[i]*sum(g*v for g, v in zip(linear, post_values[i]))
            + lam[i]*sum(g*post_values[i][r]*post_values[i][s] for r, s, g in quadratic)
            for i in range(n)
        ) % q
        assert physical_phase == low(secret)
        target = (sum(t*f(secret) for t, f in zip(shift, high_functions)) + low(secret)) % q
        assert (h(before_mask, secret) - h(mask, secret) + physical_phase - target) % q == 0
        count += 1
    return {
        "q": q, "k": k, "m": m, "ordered_points": points,
        "H_indices_zero_based": list(range(k)), "nu": nu, "lambda": lam,
        "mask_shift": shift, "linear_coefficients": linear,
        "quadratic_terms_r_s_coefficient": quadratic,
        "complete_phase_rows_checked": count,
        "checks": ["nu linear interpolation", "lambda quadratic interpolation",
                   "negative physical shift", "quadratic remainder on post-shift labels",
                   "full encoding-phase identity"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="New JSON file for detailed results.")
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error("Do not use -O: assertions perform the verification.")
    if args.output.exists():
        parser.error("Output already exists; choose a new file.")
    cases = [
        check_case(5, 2, 2, [0, 1, 3], [2, 3],
                   [lambda c: c[0]**3, lambda c: c[0]*c[1]**2],
                   [1, 2], [(0, 0, 4), (0, 1, 2), (1, 1, 3)]),
        check_case(7, 3, 1, [0, 1, 2, 4, 6], [3, 5],
                   [lambda c: c[0]**3, lambda c: c[0]**4],
                   [2], [(0, 0, 6)]),
    ]
    result = {
        "scope": "Two finite exact examples for the explicit Algorithm A formula; not a proof for arbitrary parameters.",
        "python": platform.python_version(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cases": cases,
        "total_complete_phase_rows_checked": sum(c["complete_phase_rows_checked"] for c in cases),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"cases": len(cases), "complete_phase_rows_checked": result["total_complete_phase_rows_checked"],
                      "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
