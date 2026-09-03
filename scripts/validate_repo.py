#!/usr/bin/env python3
"""Validate the manuscript-facing GMS_Ranking public release."""
from __future__ import annotations

import argparse
import hashlib
import os
import py_compile
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REQUIRED_ROOT_FILES = ["README.md", "FIGURE_MAP.tsv", "CITATION.cff"]
REQUIRED_DIRS = [
    "core", "figure1", "figure2",
    "figure3/differential_expression", "figure3/scalar_metrics", "figure3/final_assembly",
    "figure4", "figure5", "figure6_7",
    "supplement/figureS1", "supplement/figureS2", "supplement/figureS3",
    "supplement/figureS4", "supplement/figureS5",
    "config", "data", "scripts", "final_outputs",
]
FINAL_STEMS = [*(f"Figure{i}" for i in range(1, 8)), *(f"FigureS{i}" for i in range(1, 6))]
TEXT_SUFFIXES = {
    ".py", ".r", ".sh", ".md", ".txt", ".json",
    ".yaml", ".yml", ".cff", ".smk", ".gitignore",
}
SECRET_PATTERNS = [
    re.compile(r"ghp_[A-Za-z0-9]{20,}"),
    re.compile(r"github_pat_[A-Za-z0-9_]{20,}"),
    re.compile(r"AKIA[0-9A-Z]{16}"),
    re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    re.compile(r"(?i)(?:reviewer[_ ]?token|geo[_ ]?token)\s*[:=]\s*[A-Za-z0-9]{8,}"),
    re.compile(r"(?i)(?:password|passwd|api[_-]?key|secret)\s*[:=]\s*['\"]?[A-Za-z0-9_./+=-]{8,}"),
]
PRIVATE_PATTERNS = [
    re.compile(r"(?<![A-Za-z0-9_.$}])/Users/[A-Za-z0-9._-]+/"),
    re.compile(r"(?<![A-Za-z0-9_.$}])/vf/users/[A-Za-z0-9._-]+/"),
    re.compile(r"(?<![A-Za-z0-9_.$}])/home/[A-Za-z0-9._-]+/"),
    re.compile(r"(?<![A-Za-z0-9_.$}])/data/[A-Za-z0-9._-]+/(?:conda|project|work|home)/"),
]
UNFINISHED_PATTERNS = [
    re.compile(r"Implement this production script", re.I),
    re.compile(r"\bTODO\b"),
    re.compile(r"\bFIXME\b"),
    re.compile(r"MANUSCRIPT_HANDOFF", re.I),
    re.compile(r"reviewerproof", re.I),
    re.compile(r"code[-_ ]clinic", re.I),
    re.compile(r"Use only for layout/debugging", re.I),
    re.compile(r"before the public release", re.I),
]
SCAN_EXCLUSIONS = {"scripts/validate_repo.py", "SHA256SUMS.txt"}


def iter_files(root: Path) -> list[Path]:
    """Walk once without following symlinked module directories or .git."""
    files: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        dirnames[:] = [d for d in dirnames if d != ".git" and not Path(dirpath, d).is_symlink()]
        base = Path(dirpath)
        for name in filenames:
            path = base / name
            if path.is_symlink():
                continue
            files.append(path)
    return files


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--skip-checksums", action="store_true")
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    errors: list[str] = []
    warnings: list[str] = []

    if not root.is_dir():
        print(f"ERROR: repository root not found: {root}", file=sys.stderr)
        return 2

    files = iter_files(root)

    for rel in REQUIRED_ROOT_FILES:
        if not (root / rel).is_file():
            errors.append(f"missing root file: {rel}")
    for rel in REQUIRED_DIRS:
        if not (root / rel).is_dir():
            errors.append(f"missing directory: {rel}")
    for stem in FINAL_STEMS:
        if not (root / "final_outputs" / f"{stem}.pdf").is_file():
            errors.append(f"missing final display: final_outputs/{stem}.pdf")

    root_readme = (root / "README.md").resolve()
    for path in files:
        rel = path.relative_to(root).as_posix()
        size = path.stat().st_size

        if "readme" in path.name.lower() and path.resolve() != root_readme:
            errors.append(f"extra README-like file: {rel}")
        if size >= 95 * 1024 * 1024:
            errors.append(f"file >=95 MiB: {rel} ({size / 1024 / 1024:.1f} MiB)")
        elif size >= 50 * 1024 * 1024:
            warnings.append(f"large file >=50 MiB: {rel} ({size / 1024 / 1024:.1f} MiB)")

        if rel in SCAN_EXCLUSIONS or size > 5_000_000:
            continue
        scan_structured_table = path.suffix.lower() in {".tsv", ".csv"} and size <= 1_000_000
        if path.suffix.lower() in TEXT_SUFFIXES or path.name == "Makefile" or scan_structured_table:
            text = read_text(path)
            if any(pattern.search(text) for pattern in SECRET_PATTERNS):
                errors.append(f"possible credential or reviewer token: {rel}")
            if any(pattern.search(text) for pattern in PRIVATE_PATTERNS):
                errors.append(f"private absolute path: {rel}")
            if any(pattern.search(text) for pattern in UNFINISHED_PATTERNS):
                errors.append(f"unfinished or author-facing development text: {rel}")

    for forbidden in ["BUILD_PROVENANCE.txt", "PREFLIGHT_REPORT.md", "REVIEWER_GUIDE.md", "Snakefile"]:
        if (root / forbidden).exists():
            errors.append(f"forbidden root artifact: {forbidden}")
    for rel in [
        "supplement/figureS1/scripts/run.py",
        "supplement/figureS4/scripts/run.py",
        "figure2/reviewerproof",
    ]:
        if (root / rel).exists():
            errors.append(f"unfinished/development path retained: {rel}")

    expected_links = {
        "supplement/figureS2/code": "../../figure5",
        "supplement/figureS3/code": "../../figure6_7",
        "supplement/figureS5/code": "../../figure3/differential_expression",
    }
    for rel, target in expected_links.items():
        path = root / rel
        if not path.is_symlink():
            errors.append(f"missing symlink: {rel}")
        elif path.readlink().as_posix() != target or not path.exists():
            errors.append(f"bad symlink: {rel} -> {path.readlink()}")

    trace = root / "supplement/figureS1/source_data/GEO_TRACE.tsv"
    trace_text = read_text(trace)
    if not all(token in trace_text for token in ["GSE318271", "GSE318584", "FigureS1A", "FigureS1B"]):
        errors.append("invalid Figure S1 GEO trace")

    # Compile all Python and Bash entry points.
    for path in files:
        rel = path.relative_to(root).as_posix()
        if path.suffix.lower() == ".py":
            try:
                py_compile.compile(str(path), doraise=True)
            except Exception as exc:
                errors.append(f"Python syntax failure: {rel}: {exc}")
        elif path.suffix.lower() == ".sh":
            result = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
            if result.returncode:
                errors.append(f"Bash syntax failure: {rel}: {result.stderr.strip()}")

    # Parse all R files with a temporary helper script when Rscript is available.
    # Avoid `Rscript -e ... <many files>` because some macOS/R installations
    # misparse the expression and report spurious `ARGUMENT ... __ignored__` errors.
    r_files = [str(p) for p in files if p.suffix.lower() == ".r"]
    if r_files and shutil.which("Rscript"):
        with tempfile.TemporaryDirectory(prefix="gms_r_syntax_") as tmp_dir:
            tmp = Path(tmp_dir)
            file_list = tmp / "r_files.txt"
            helper = tmp / "parse_all.R"
            file_list.write_text("\n".join(r_files) + "\n", encoding="utf-8")
            helper.write_text(
                "args <- commandArgs(trailingOnly=TRUE)\n"
                "if (length(args) != 1L) stop('Expected one R-file-list argument')\n"
                "files <- readLines(args[[1]], warn=FALSE)\n"
                "failed <- FALSE\n"
                "for (f in files) {\n"
                "  tryCatch(parse(file=f), error=function(e) {\n"
                "    cat(sprintf('PARSE_FAIL\\t%s\\t%s\\n', f, conditionMessage(e)))\n"
                "    failed <<- TRUE\n"
                "  })\n"
                "}\n"
                "quit(status=if (failed) 1L else 0L)\n",
                encoding="utf-8",
            )
            result = subprocess.run(
                ["Rscript", "--vanilla", str(helper), str(file_list)],
                capture_output=True,
                text=True,
            )
        if result.returncode:
            errors.append("R syntax validation failed: " + (result.stdout + result.stderr).strip())

    expected_header = "final_display\tpanels\tanalysis_directory\tprimary_code\tprimary_source_data\tfinal_output"
    map_path = root / "FIGURE_MAP.tsv"
    if map_path.is_file():
        lines = map_path.read_text(encoding="utf-8", errors="replace").splitlines()
        if not lines or lines[0] != expected_header:
            errors.append("FIGURE_MAP.tsv header mismatch")

    # No two manuscript display PDFs should be byte-identical.
    hashes: dict[str, list[str]] = {}
    for stem in FINAL_STEMS:
        path = root / "final_outputs" / f"{stem}.pdf"
        if path.is_file():
            hashes.setdefault(sha256(path), []).append(path.name)
    for names in hashes.values():
        if len(names) > 1:
            errors.append("duplicate final PDF content: " + ", ".join(sorted(names)))

    if not args.skip_checksums:
        checksum_lines = []
        for path in sorted(files):
            rel = path.relative_to(root).as_posix()
            if rel != "SHA256SUMS.txt":
                checksum_lines.append(f"{sha256(path)}  {rel}")
        (root / "SHA256SUMS.txt").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")

    print(f"Validation complete: {len(errors)} error(s), {len(warnings)} warning(s)")
    for item in errors:
        print("ERROR: " + item, file=sys.stderr)
    for item in warnings:
        print("WARNING: " + item, file=sys.stderr)
    return 1 if errors or (args.strict and warnings) else 0


if __name__ == "__main__":
    raise SystemExit(main())

