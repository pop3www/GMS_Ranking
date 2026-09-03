from __future__ import annotations

import json
import platform
import sys


def main() -> int:
    try:
        import matplotlib
        import numpy
        import pandas
        import scipy
        import sklearn
    except Exception as exc:  # pragma: no cover - diagnostic path
        sys.stderr.write("[ERROR] Figure 6 environment preflight failed while importing core packages.\n")
        sys.stderr.write(f"[ERROR] {type(exc).__name__}: {exc}\n")
        sys.stderr.write(
            "[ERROR] This usually means the active Python environment is incompatible with the pinned bundle "
            "versions (for example numpy / scikit-learn binary mismatch). Activate fig6_repro or another clean env "
            "built from environment.yml, then rerun.\n"
        )
        return 1

    payload = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "python_executable": sys.executable,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    print("[OK] Figure 6 environment preflight passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
