#!/usr/bin/env python3
"""Compatibility wrapper for the Fig. 7/8 R workflow.

The repository scaffold may point at scripts/run.py.  The production workflow is
implemented in scripts/run.R, so this wrapper simply delegates to Rscript while
preserving arguments.
"""
from __future__ import annotations
import os
import subprocess
import sys
from pathlib import Path


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    run_r = script_dir / "run.R"
    if not run_r.exists():
        sys.stderr.write(f"ERROR: expected {run_r}\n")
        return 2
    cmd = [os.environ.get("RSCRIPT", "Rscript"), str(run_r), *sys.argv[1:]]
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
