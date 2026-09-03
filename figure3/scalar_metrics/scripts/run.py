#!/usr/bin/env python3
"""Python wrapper for the Figure 4 R production script."""
from __future__ import annotations
import pathlib
import subprocess
import sys

script_dir = pathlib.Path(__file__).resolve().parent
r_script = script_dir / "run.R"
if not r_script.exists():
    raise SystemExit(f"Missing R script: {r_script}")
cmd = ["Rscript", str(r_script), *sys.argv[1:]]
raise SystemExit(subprocess.call(cmd))
