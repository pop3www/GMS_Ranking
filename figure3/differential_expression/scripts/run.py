#!/usr/bin/env python3
"""Python wrapper for the Figure 3 R production entry point."""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

def main() -> int:
    script = Path(__file__).with_name("run.R")
    if not script.exists():
        sys.stderr.write(f"Missing R script: {script}\n")
        return 1
    cmd = ["Rscript", str(script), *sys.argv[1:]]
    return subprocess.call(cmd)

if __name__ == "__main__":
    raise SystemExit(main())
