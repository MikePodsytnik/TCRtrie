import os
import subprocess
import sys
from pathlib import Path

def _binary_path() -> Path:
    package_dir = Path(__file__).resolve().parent
    candidates = [
        package_dir / "tcrtrie",
        package_dir / "tcrtrie.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise RuntimeError("Bundled tcrtrie executable was not found")

def main() -> None:
    binary = _binary_path()
    completed = subprocess.run([str(binary), *sys.argv[1:]])
    raise SystemExit(completed.returncode)

if __name__ == "__main__":
    main()
