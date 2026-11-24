#!/usr/bin/env python3
"""Run the fixed finite checks and compare their mathematical output to baselines."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
ABS_TOL = 1e-12
REL_TOL = 1e-10


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path):
    def reject_constant(value):
        raise ValueError(f"Non-finite JSON number: {value}")
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject_constant)


def mathematical_result(data, excluded_fields):
    """Remove only explicitly listed top-level runtime/provenance fields."""
    if isinstance(data, dict):
        data = {key: value for key, value in data.items() if key not in excluded_fields}
        # LAPACK may round the same zero eigenvalue to '-0.0' or '0.0'.
        if "equal_weight_spectrum" in data:
            spectrum = {}
            for key, multiplicity in data["equal_weight_spectrum"].items():
                number = float(key)
                if not math.isfinite(number):
                    raise ValueError("Non-finite spectrum key")
                key = str(0.0 if number == 0 else number)
                spectrum[key] = spectrum.get(key, 0) + multiplicity
            data["equal_weight_spectrum"] = spectrum
    return data


def differences(actual, expected, path="$", found=None):
    """Compare every key and item; integers are exact, floats use stated tolerances."""
    found = [] if found is None else found
    if isinstance(expected, bool):
        if type(actual) is not bool or actual != expected:
            found.append(f"{path}: expected {expected!r}, got {actual!r}")
    elif isinstance(expected, int):
        if type(actual) is not int or actual != expected:
            found.append(f"{path}: expected integer {expected}, got {actual!r}")
    elif isinstance(expected, float):
        if (type(actual) not in (int, float) or not math.isfinite(actual)
                or not math.isclose(actual, expected, rel_tol=REL_TOL, abs_tol=ABS_TOL)):
            found.append(f"{path}: expected {expected!r}, got {actual!r}")
    elif isinstance(expected, dict):
        if not isinstance(actual, dict):
            found.append(f"{path}: expected object")
        else:
            if actual.keys() != expected.keys():
                found.append(f"{path}: key mismatch; missing={sorted(expected.keys()-actual.keys())}, "
                             f"extra={sorted(actual.keys()-expected.keys())}")
            for key in expected.keys() & actual.keys():
                differences(actual[key], expected[key], f"{path}.{key}", found)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            found.append(f"{path}: expected list of length {len(expected)}")
        else:
            for i, (left, right) in enumerate(zip(actual, expected)):
                differences(left, right, f"{path}[{i}]", found)
    elif type(actual) is not type(expected) or actual != expected:
        found.append(f"{path}: expected {expected!r}, got {actual!r}")
    return found


def validate_diagnostics(check_id, data):
    """Explicit acceptance of diagnostics whose original script only printed values."""
    errors = []
    def require(condition, message):
        if not condition:
            errors.append(message)
    def near(value, target=0):
        return math.isfinite(value) and abs(value - target) <= ABS_TOL
    if check_id == "E03":
        require(data["noncanonical_choi_rank"] == 2, "E03 Choi rank != 2")
        for key in ("choi_reference_marginal_error", "unrecorded_twirl_magic_fixed_error",
                    "unrecorded_twirl_offdiagonal_norm", "key_view_max_difference",
                    "cipher_view_max_difference"):
            require(near(data[key]), f"E03 {key} exceeds tolerance")
        require(all(near(x, 1) for x in data["all_secret_recovery_probabilities"]),
                "E03 recovery probabilities != 1")
    elif check_id == "E04":
        require(data["choi_rank"] == 2, "E04 Choi rank != 2")
        spectrum = sorted(data["choi_eigenvalues"])
        require(len(spectrum) == 25 and all(near(x) for x in spectrum[:23])
                and all(near(x, 0.5) for x in spectrum[23:]), "E04 Choi spectrum mismatch")
        require(near(data["choi_trace"], 1), "E04 trace != 1")
        require(near(data["output_wigner_sum"], 1), "E04 Wigner sum != 1")
        require(near(data["output_wigner_min"], -math.sqrt(5)/25), "E04 Wigner minimum mismatch")
        require(near(data["max_branch_effect_error"]), "E04 branch effect error exceeds tolerance")
    elif check_id == "E06":
        require(data["polynomial_direction_ranks"] == [2, 3, 4], "E06 direction ranks mismatch")
        require(data["logical_counts_from_direction_ranks"] == [25, 125, 625], "E06 library counts mismatch")
        for key in ("D27_yB_zero_mod5", "D27_yJ_zero_mod5", "light_shift_identity_all_zero",
                    "D10_E5_match", "D11_U012_match", "D12_T34_match"):
            require(data[key] is True, f"E06 {key} must be true")
        # This is deliberately FALSE: the modular certificate is not an integer annihilator.
        require(data["D27_yB_exact_integer_zero"] is False,
                "E06 modular and integer certificates were conflated")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="New directory; existing paths are refused.")
    parser.add_argument("--checks", nargs="+", help="Selected check IDs (default: all ten).")
    parser.add_argument("--list", action="store_true", help="List check IDs without executing them.")
    parser.add_argument("--timeout", type=float, default=180, help="Seconds per script (default: 180).")
    args = parser.parse_args()
    specs = load_json(ROOT / "inputs/checks.json")
    if args.list:
        for spec in specs:
            print(f"{spec['id']}: {spec['description']}")
        return 0
    if sys.flags.optimize:
        parser.error("Do not use -O or PYTHONOPTIMIZE: the script assertions are part of verification.")
    if sys.version_info < (3, 12):
        parser.error("This package requires Python 3.12 or newer; tested with Python 3.12.14.")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be a positive finite number")
    selected = set(args.checks or [spec["id"] for spec in specs])
    unknown = selected - {spec["id"] for spec in specs}
    if unknown:
        parser.error(f"Unknown checks: {sorted(unknown)}")
    specs = [spec for spec in specs if spec["id"] in selected]
    try:
        numpy_version = importlib.metadata.version("numpy")
    except importlib.metadata.PackageNotFoundError:
        parser.error("NumPy is missing. Install dependencies with: python -m pip install -r requirements.txt")
    output = (args.output_dir or ROOT / "results" / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")).resolve()
    if output == ROOT or ROOT.is_relative_to(output):
        parser.error("Output must not be the repository or its parent.")
    # Refuse reuse even when an earlier run has failed: every run gets its own evidence.
    if output.exists():
        parser.error(f"Output already exists; select a new directory: {output}")
    output.mkdir(parents=True)
    manifest = {
        "package_version": (ROOT / "VERSION").read_text().strip(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "numpy": numpy_version,
        "platform": platform.system(), "machine": platform.machine(),
        "float_absolute_tolerance": ABS_TOL, "float_relative_tolerance": REL_TOL,
        "exact_integer_comparison": True, "assertions_enabled": True,
        "reference_policy": "Archived mathematical results; runtime/provenance exclusions are explicit in inputs/checks.json.",
        "checks": [],
    }
    environment = os.environ.copy()
    for key in ("PYTHONOPTIMIZE", "PYTHONPATH", "THESIS_WEB_AUDIT_SOURCE", "THESIS_AUDIT_SOURCE_ROOT",
                "THESIS_AUDIT_OUTPUT_DIR"):
        environment.pop(key, None)
    environment.update(THESIS_WEB_AUDIT_RESULTS=".", PYTHONHASHSEED="0",
                       OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for spec in specs:
        case_dir = output / spec["id"]
        case_dir.mkdir()
        command = [sys.executable, "-I", "-B", str(ROOT / spec["script"]), *spec.get("args", [])]
        entry = {"id": spec["id"], "script": spec["script"],
                 "script_sha256": sha256(ROOT / spec["script"]),
                 "reference_sha256": sha256(ROOT / spec["reference"])}
        start = time.monotonic()
        try:
            completed = subprocess.run(command, cwd=case_dir, env=environment, capture_output=True,
                                       text=True, encoding="utf-8", timeout=args.timeout)
            (case_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
            (case_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
            entry["returncode"] = completed.returncode
            if completed.returncode != 0:
                raise RuntimeError(f"Script exited with code {completed.returncode}; see stderr.txt")
            raw_path = case_dir / "raw.json"
            if spec["result"] == "stdout":
                raw_path.write_text(completed.stdout, encoding="utf-8")
            else:
                raw_path = case_dir / spec["result"]
            data = mathematical_result(load_json(raw_path), spec.get("exclude_fields", []))
            expected = load_json(ROOT / spec["reference"])
            errors = differences(data, expected) + validate_diagnostics(spec["id"], data)
            (case_dir / "mathematical_result.json").write_text(
                json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            entry.update(status="PASS" if not errors else "FAIL", mismatches=errors,
                         result_sha256=sha256(case_dir / "mathematical_result.json"))
        except Exception as exc:
            entry.update(status="FAIL", error=f"{type(exc).__name__}: {exc}")
        entry["elapsed_seconds"] = round(time.monotonic() - start, 4)
        manifest["checks"].append(entry)
        print(f"{entry['status']} {spec['id']} ({entry['elapsed_seconds']:.3f}s)", flush=True)
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["passed"] = sum(x["status"] == "PASS" for x in manifest["checks"])
    manifest["failed"] = len(specs) - manifest["passed"]
    (output / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Results: {output}")
    return 0 if manifest["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
