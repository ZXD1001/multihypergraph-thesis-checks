#!/usr/bin/env python3
"""Verify the frozen source, documentation, and baseline files without dependencies."""
import hashlib
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parent
manifest = json.loads((root / "SHA256SUMS.json").read_text(encoding="utf-8"))
failures = []
for name, expected in manifest["files"].items():
    path = root / name
    if path.is_symlink() or not path.is_file():
        failures.append(f"Missing/nonregular: {name}")
    elif hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        failures.append(f"Changed: {name}")
if failures:
    print("\n".join(failures))
    sys.exit(1)
print(f"PASS: {len(manifest['files'])} frozen files match version {manifest['version']}")
