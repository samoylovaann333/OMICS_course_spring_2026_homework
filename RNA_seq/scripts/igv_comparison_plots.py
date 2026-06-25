#!/usr/bin/env python3
"""Render the bigWig comparisons asked for in the IGV task as static figures
(stand-in for IGV screenshots, since pyBigWig + matplotlib reproduce the same
coverage tracks programmatically and the plots are easy to regenerate)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyBigWig

ROOT = Path(__file__).resolve().parent.parent
TRACKS = ROOT / "results" / "tracks"
FIGURES = ROOT / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

SAMPLES = ["MoPh7", "MoPh11", "MoPh14", "MoPh15"]
COLORS = {"MoPh7": "#1f77b4", "MoPh11": "#ff7f0e", "MoPh14": "#2ca02c", "MoPh15": "#d62728"}

DIFF_REGION = ("chr1", 35_491_500, 35_495_500)
DIFF_BIN = ("chr1", 35_492_300, 35_494_300)
ZOOM_REGION = ("chr1", 35_494_280, 35_494_360)


def track_values(path: Path, chrom: str, start: int, end: int, nbins: int) -> np.ndarray:
    bw = pyBigWig.open(str(path))
    vals = bw.stats(chrom, start, end, nBins=nbins, type="mean")
    bw.close()
    return np.nan_to_num(np.array(vals, dtype=float))


def plot_rnaseq_vs_cage_across_samples() -> None:
    chrom, start, end = DIFF_REGION
    nbins = 400
    x = np.linspace(start, end, nbins)

    fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    for sample in SAMPLES:
        rnaseq = track_values(TRACKS / "rnaseq" / f"{sample}.rnaseq.STAR.bw", chrom, start, end, nbins)
        axes[0].plot(x, rnaseq, label=sample, color=COLORS[sample], lw=1.2)
        cage = track_values(TRACKS / "cage" / f"{sample}.cage.STAR.bw", chrom, start, end, nbins)
        axes[1].plot(x, cage, label=sample, color=COLORS[sample], lw=1.2)

    axes[0].set_ylabel("RNA-seq coverage")
    axes[1].set_ylabel("CAGE coverage")
    axes[1].set_xlabel(f"{chrom} position")
    axes[0].axvspan(DIFF_BIN[1], DIFF_BIN[2], color="grey", alpha=0.2)
    axes[1].axvspan(DIFF_BIN[1], DIFF_BIN[2], color="grey", alpha=0.2)
    axes[0].set_title(f"RNA-seq vs CAGE, {chrom}:{start}-{end} (4 samples, STAR alignments)")
    axes[0].legend(ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "rnaseq_vs_cage_chr15_region.png", dpi=150)
    plt.close(fig)


def plot_star_vs_bwa_moph7() -> None:
    chrom, start, end = ZOOM_REGION
    nbins = end - start  # 1 bp resolution to resolve the junction ramp

    star = track_values(TRACKS / "rnaseq" / "MoPh7.rnaseq.STAR.bw", chrom, start, end, nbins)
    bwa = track_values(TRACKS / "bwa" / "MoPh7.rnaseq.BWA.q30.bw", chrom, start, end, nbins)
    x = np.arange(start, end)

    fig, ax = plt.subplots(figsize=(11, 4))
    ax.step(x, star, where="mid", label="STAR (splice-aware)", color="#1f77b4", lw=1.8)
    ax.step(x, bwa, where="mid", label="BWA (no splice junctions)", color="#d62728", lw=1.8)
    ax.set_title(f"MoPh7: STAR vs BWA at an exon boundary, {chrom}:{start}-{end}")
    ax.set_xlabel(f"{chrom} position")
    ax.set_ylabel("Coverage")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "star_vs_bwa_moph7_chr15_zoom.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    plot_rnaseq_vs_cage_across_samples()
    plot_star_vs_bwa_moph7()
    print("Saved figures to", FIGURES)
