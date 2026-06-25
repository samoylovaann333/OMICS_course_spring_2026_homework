#!/usr/bin/env python3
"""Scan RNA-seq STAR tracks of all 4 samples for a genomic region where
signal differs notably between samples (used to pick a region for the IGV task)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pyBigWig

SAMPLES = ["MoPh7", "MoPh11", "MoPh14", "MoPh15"]
BIN_SIZE = 5_000
TRACKS_DIR = Path(__file__).resolve().parent.parent / "results" / "tracks" / "rnaseq"


def main() -> None:
    bws = {s: pyBigWig.open(str(TRACKS_DIR / f"{s}.rnaseq.STAR.bw")) for s in SAMPLES}
    chroms = bws["MoPh7"].chroms()

    candidates = []  # (score, chrom, start, end, means)
    for chrom, length in chroms.items():
        nbins = length // BIN_SIZE
        if nbins == 0:
            continue
        per_sample = []
        for s in SAMPLES:
            vals = bws[s].stats(chrom, 0, nbins * BIN_SIZE, nBins=nbins, type="mean")
            per_sample.append(np.nan_to_num(np.array(vals, dtype=float)))
        mat = np.vstack(per_sample)  # 4 x nbins

        mean_signal = mat.mean(axis=0)
        std_signal = mat.std(axis=0)
        cv = np.divide(std_signal, mean_signal, out=np.zeros_like(std_signal), where=mean_signal > 0)
        # moderate absolute signal (avoid both noise floor and single-bp mega-spikes)
        mask = (mean_signal > 5) & (mean_signal < 2000)
        if not mask.any():
            continue
        score = np.where(mask, cv * np.sqrt(mean_signal), 0)
        idx = int(np.argmax(score))
        candidates.append((score[idx], chrom, idx * BIN_SIZE, (idx + 1) * BIN_SIZE, mat[:, idx]))

    for bw in bws.values():
        bw.close()

    candidates.sort(key=lambda x: -x[0])
    print("Top candidate regions (score = CV * sqrt(mean), moderate-signal bins):\n")
    for score, chrom, start, end, means in candidates[:10]:
        win_start, win_end = max(0, start - 20_000), end + 20_000
        means_str = ", ".join(f"{s}={v:.1f}" for s, v in zip(SAMPLES, means))
        print(f"{chrom}:{win_start}-{win_end}  score={score:.2f}  [{means_str}]")


if __name__ == "__main__":
    main()
