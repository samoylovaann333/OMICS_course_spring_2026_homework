#!/usr/bin/env python3
"""Save insulation-score domain boundaries as bedGraph (insulation_and_boundaries.ipynb,
'Задание для самостоятельного разбора'), for later intersection with CTCF on the ChIP-seq day.
"""

from __future__ import annotations

from pathlib import Path

import cooler
import pandas as pd
from cooltools import insulation

ROOT = Path(__file__).resolve().parent.parent
MCOOL_PATH = "/Users/anna/Desktop/day2_HiC_practice/data/MoPh7_enr_v2.mcool"
RESOLUTION = 100_000
WINDOW = 300_000
SAMPLE = "MoPh7_enr_v2"


def main() -> None:
    clr = cooler.Cooler(f"{MCOOL_PATH}::/resolutions/{RESOLUTION}")

    canonical_chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    chromsizes = clr.chromsizes.reset_index()
    chromsizes.columns = ["chrom", "length"]
    chromsizes = chromsizes[chromsizes["chrom"].isin(canonical_chroms)].copy()

    view_df = pd.DataFrame(
        {
            "chrom": chromsizes["chrom"],
            "start": 0,
            "end": chromsizes["length"],
            "name": chromsizes["chrom"],
        }
    )

    insulation_table = insulation(
        clr, [WINDOW], view_df=view_df, clr_weight_name="weight", nproc=4
    )

    boundary_col = f"is_boundary_{WINDOW}"
    strength_col = f"boundary_strength_{WINDOW}"

    boundaries = insulation_table[insulation_table[boundary_col] == True][  # noqa: E712
        ["chrom", "start", "end", strength_col]
    ]

    out_dir = ROOT / "results" / "boundaries"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{SAMPLE}_boundaries_{WINDOW // 1000}kb.bedGraph"
    boundaries.to_csv(out_path, sep="\t", header=False, index=False)
    print(f"Saved {len(boundaries)} boundaries to {out_path}")


if __name__ == "__main__":
    main()
