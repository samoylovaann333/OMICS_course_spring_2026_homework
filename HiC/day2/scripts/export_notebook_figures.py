#!/usr/bin/env python3
"""Dump embedded image outputs from the day2 notebooks into results/figures/
so they can be referenced in README/slides instead of staying buried in .ipynb."""

from __future__ import annotations

import base64
import json
from pathlib import Path

NOTEBOOKS_DIR = Path(__file__).resolve().parent.parent / "notebooks"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "results" / "figures"


def export(notebook: Path) -> None:
    nb = json.loads(notebook.read_text())
    stem = notebook.stem
    img_i = 0
    for cell in nb["cells"]:
        for output in cell.get("outputs", []):
            data = output.get("data", {})
            for mime, payload in data.items():
                if not mime.startswith("image/"):
                    continue
                ext = mime.split("/")[-1]
                img_i += 1
                out_path = FIGURES_DIR / f"{stem}_{img_i:02d}.{ext}"
                if isinstance(payload, list):
                    payload = "".join(payload)
                out_path.write_bytes(base64.b64decode(payload))
                print("wrote", out_path)


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for nb_path in sorted(NOTEBOOKS_DIR.glob("*.ipynb")):
        export(nb_path)


if __name__ == "__main__":
    main()
