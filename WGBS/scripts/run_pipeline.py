#!/usr/bin/env python3
"""
День 4. WGBS/Bismark: анализ метилирования
==========================================
Запуск: python3 run_pipeline.py

Скрипт последовательно выполняет все три этапа анализа:
  1. QC и подсчёт beta-value / M-value
  2. Создание bigWig треков для IGV
  3. Интеграция метилирования с ChIP-seq (если есть файлы)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid")

# =============================================================================
# ПУТИ
# =============================================================================

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "data" / "bismark"
TABLES_DIR = ROOT / "results" / "tables"
TRACKS_DIR = ROOT / "results" / "tracks"
FIGURES_DIR = ROOT / "results" / "figures"

REF_DIR = ROOT.parent / "day1_HiC_practice" / "data" / "reference"
CHROM_SIZES_FILE = REF_DIR / "chrom.sizes"
REFERENCE_FASTA = REF_DIR / "T2T_human.fna"

BISMARK_COV = DATA_DIR / "MoPh7_1_bismark_bt2_pe.bismark.cov.gz"
METHYLATION_TABLE = TABLES_DIR / "MoPh7_cpg_methylation_values.tsv.gz"

CHIP_TRACKS = {
    "H3K27ac": ROOT.parent / "day3_ChIPseq_practice" / "results" / "macs" / "MoPh7_H3K27Ac" / "MoPh7_H3K27Ac_FE.bw",
    "H3K9me3": ROOT.parent / "day3_ChIPseq_practice" / "results" / "macs" / "MoPh7_H3K9me3" / "MoPh7_H3K9me3_FE.bw",
}

SAMPLE = "MoPh7"
MIN_COVERAGE = 5      # минимальное покрытие CpG
MAX_ROWS_QC = 1_000_000  # строк для QC-графиков (весь файл — None)

for d in [TABLES_DIR, TRACKS_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# =============================================================================
# ЧАСТЬ 1 — QC и подсчёт метилирования
# =============================================================================

def run_part1():
    print("\n" + "=" * 60)
    print("ЧАСТЬ 1: QC и подсчёт beta-value / M-value")
    print("=" * 60)

    if not BISMARK_COV.exists():
        print(f"ERROR: файл не найден: {BISMARK_COV}")
        print("Запусти сначала: python3 scripts/download_bismark_data.py")
        sys.exit(1)

    print(f"Читаем {BISMARK_COV.name} (первые {MAX_ROWS_QC:,} строк для QC)...")

    df_qc = pd.read_csv(
        BISMARK_COV,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "meth_percent", "meth_count", "unmeth_count"],
        nrows=MAX_ROWS_QC,
    )

    df_qc["coverage"] = df_qc["meth_count"] + df_qc["unmeth_count"]
    df_qc["beta_value"] = df_qc["meth_count"] / df_qc["coverage"]
    df_qc["m_value"] = np.log2(
        (df_qc["meth_count"] + 1) / (df_qc["unmeth_count"] + 1)
    )

    print(f"  Всего CpG в подвыборке:  {len(df_qc):,}")
    print(f"  Медианное покрытие:       {df_qc['coverage'].median():.1f}")
    print(f"  CpG с coverage >= {MIN_COVERAGE}:    {(df_qc['coverage'] >= MIN_COVERAGE).sum():,}")

    # -- График 1: распределение покрытия --
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].hist(df_qc["coverage"].clip(upper=50), bins=50, color="steelblue", edgecolor="white")
    axes[0].axvline(MIN_COVERAGE, color="red", linestyle="--", label=f"min coverage = {MIN_COVERAGE}")
    axes[0].set_xlabel("Coverage")
    axes[0].set_ylabel("Number of CpGs")
    axes[0].set_title("Coverage distribution")
    axes[0].legend()

    df_filtered_qc = df_qc[df_qc["coverage"] >= MIN_COVERAGE]

    axes[1].hist(df_filtered_qc["beta_value"], bins=80, color="mediumseagreen", edgecolor="white")
    axes[1].set_xlabel("Beta-value")
    axes[1].set_ylabel("Number of CpGs")
    axes[1].set_title(f"{SAMPLE}: Beta-value distribution\n(coverage ≥ {MIN_COVERAGE})")

    axes[2].hist(df_filtered_qc["m_value"], bins=80, color="coral", edgecolor="white")
    axes[2].set_xlabel("M-value")
    axes[2].set_ylabel("Number of CpGs")
    axes[2].set_title(f"{SAMPLE}: M-value distribution\n(coverage ≥ {MIN_COVERAGE})")

    plt.tight_layout()
    fig_path = FIGURES_DIR / f"{SAMPLE}_methylation_distributions.png"
    plt.savefig(fig_path, dpi=150)
    plt.close()
    print(f"  Сохранён график: {fig_path}")

    # -- Читаем весь файл и сохраняем таблицу --
    print("Читаем весь файл и фильтруем...")
    df_full = pd.read_csv(
        BISMARK_COV,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "meth_percent", "meth_count", "unmeth_count"],
    )
    df_full["coverage"] = df_full["meth_count"] + df_full["unmeth_count"]
    df_full["beta_value"] = df_full["meth_count"] / df_full["coverage"]
    df_full["m_value"] = np.log2(
        (df_full["meth_count"] + 1) / (df_full["unmeth_count"] + 1)
    )

    filtered = df_full[df_full["coverage"] >= MIN_COVERAGE].copy()
    print(f"  После фильтрации: {len(filtered):,} CpG")

    filtered.to_csv(METHYLATION_TABLE, sep="\t", index=False, compression="gzip")
    print(f"  Таблица сохранена: {METHYLATION_TABLE}")

    return filtered


# =============================================================================
# ЧАСТЬ 2 — bigWig треки
# =============================================================================

def read_chrom_sizes(path):
    chrom_sizes = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def write_bigwig(table, value_col, output_path, chrom_sizes):
    import pyBigWig

    t = table[["chrom", "bw_start", "bw_end", value_col]].copy()
    t = t.replace([np.inf, -np.inf], np.nan).dropna()
    t = t[t["chrom"].isin(chrom_sizes)]

    # Сортируем по порядку хромосом как в chrom_sizes
    chrom_order = {c: i for i, c in enumerate(chrom_sizes)}
    t["_order"] = t["chrom"].map(chrom_order)
    t = t.sort_values(["_order", "bw_start"]).drop(columns="_order")

    bw = pyBigWig.open(str(output_path), "w")
    bw.addHeader(list(chrom_sizes.items()))

    for chrom, group in t.groupby("chrom", sort=False):
        bw.addEntries(
            [chrom] * len(group),
            group["bw_start"].tolist(),
            ends=group["bw_end"].tolist(),
            values=group[value_col].astype(float).tolist(),
        )

    bw.close()
    return output_path


def run_part2(filtered):
    print("\n" + "=" * 60)
    print("ЧАСТЬ 2: Создание bigWig треков")
    print("=" * 60)

    if not CHROM_SIZES_FILE.exists():
        print(f"WARNING: chrom.sizes не найден: {CHROM_SIZES_FILE}")
        print("Пропускаем часть 2. Запусти шаг 6 из инструкции.")
        return

    chrom_sizes = read_chrom_sizes(CHROM_SIZES_FILE)
    print(f"  Загружено хромосом: {len(chrom_sizes)}")

    df = filtered.copy()
    df["bw_start"] = (df["start"] - 1).clip(lower=0).astype(int)
    df["bw_end"] = df["bw_start"] + 1
    df = df[df["chrom"].isin(chrom_sizes)]
    df = (
        df.groupby(["chrom", "bw_start", "bw_end"], as_index=False)
        .agg({"beta_value": "mean", "m_value": "mean", "coverage": "mean"})
    )

    tracks = {
        "beta_methylation": "beta_value",
        "m_value": "m_value",
        "coverage": "coverage",
    }

    for track_name, col in tracks.items():
        bedgraph_path = TRACKS_DIR / f"{SAMPLE}_{track_name}.bedGraph"
        bigwig_path = TRACKS_DIR / f"{SAMPLE}_{track_name}.bw"

        bed = df[["chrom", "bw_start", "bw_end", col]].copy()
        bed = bed.replace([np.inf, -np.inf], np.nan).dropna()
        bed.to_csv(bedgraph_path, sep="\t", header=False, index=False)

        write_bigwig(df, col, bigwig_path, chrom_sizes)
        print(f"  Создан трек: {bigwig_path.name}")

    # GC content и CpG obs/exp (нужен референс)
    if not REFERENCE_FASTA.exists():
        print(f"  WARNING: референс не найден, пропускаем GC/CpG треки")
        return

    print("  Считаем GC content и CpG obs/exp (это займёт несколько минут)...")
    import pysam

    bin_size = 100
    fasta = pysam.FastaFile(str(REFERENCE_FASTA))

    regions = (
        df.groupby("chrom")
        .agg(start=("bw_start", "min"), end=("bw_end", "max"))
        .reset_index()
    )

    gc_rows = []
    for _, row in regions.iterrows():
        chrom = row["chrom"]
        chrom_size = chrom_sizes.get(chrom, 0)
        region_start = int(row["start"] // bin_size * bin_size)
        region_end = int(np.ceil(row["end"] / bin_size) * bin_size)
        region_end = min(region_end, chrom_size)

        for start in range(region_start, region_end, bin_size):
            end = min(start + bin_size, chrom_size)
            seq = fasta.fetch(chrom, start, end).upper()
            length = len(seq)
            if length == 0:
                continue

            c_count = seq.count("C")
            g_count = seq.count("G")
            cg_count = seq.count("CG")
            gc_content = (c_count + g_count) / length
            expected_cpg = (c_count * g_count) / length if length else 0
            cpg_oe = cg_count / expected_cpg if expected_cpg > 0 else np.nan

            gc_rows.append((chrom, start, end, gc_content, cpg_oe))

    gc_df = pd.DataFrame(
        gc_rows, columns=["chrom", "bw_start", "bw_end", "gc_content", "cpg_obs_exp"]
    )

    for track_name, col in {"gc_content": "gc_content", "cpg_obs_exp": "cpg_obs_exp"}.items():
        bedgraph_path = TRACKS_DIR / f"T2T_{track_name}_100bp.bedGraph"
        bigwig_path = TRACKS_DIR / f"T2T_{track_name}_100bp.bw"

        bed = gc_df[["chrom", "bw_start", "bw_end", col]].copy()
        bed = bed.replace([np.inf, -np.inf], np.nan).dropna()
        bed.to_csv(bedgraph_path, sep="\t", header=False, index=False)

        write_bigwig(gc_df, col, bigwig_path, chrom_sizes)
        print(f"  Создан трек: {bigwig_path.name}")


# =============================================================================
# ЧАСТЬ 3 — интеграция с ChIP-seq
# =============================================================================

def run_part3():
    print("\n" + "=" * 60)
    print("ЧАСТЬ 3: Интеграция с ChIP-seq")
    print("=" * 60)

    import pyBigWig
    from scipy.stats import pearsonr, spearmanr

    missing = [name for name, path in CHIP_TRACKS.items() if not path.exists()]
    if missing:
        print(f"  WARNING: не найдены ChIP-seq треки: {missing}")
        print("  Пропускаем часть 3.")
        return

    print("  Читаем таблицу метилирования...")
    df = pd.read_csv(METHYLATION_TABLE, sep="\t")

    max_rows = 1_000_000
    if len(df) > max_rows:
        df = df.sample(max_rows, random_state=1).sort_values(["chrom", "start", "end"])

    print(f"  CpG для анализа: {len(df):,}")

    def mean_bigwig_signal(bw, chrom, start, end, window=500):
        if chrom not in bw.chroms():
            return None
        chrom_size = bw.chroms()[chrom]
        center = int((start + end) / 2)
        left = max(0, center - window)
        right = min(chrom_size, center + window)
        return bw.stats(chrom, left, right, type="mean")[0]

    for mark_name, bw_path in CHIP_TRACKS.items():
        print(f"  Извлекаем сигнал {mark_name}...")
        bw = pyBigWig.open(str(bw_path))
        df[mark_name] = [
            mean_bigwig_signal(bw, row.chrom, row.start, row.end)
            for row in df.itertuples()
        ]
        bw.close()

    df_clean = df.dropna(subset=list(CHIP_TRACKS.keys()) + ["beta_value"])
    print(f"  CpG с сигналом для корреляции: {len(df_clean):,}")

    # Корреляции
    results = []
    for mark_name in CHIP_TRACKS:
        r_p, p_p = pearsonr(df_clean["beta_value"], df_clean[mark_name])
        r_s, p_s = spearmanr(df_clean["beta_value"], df_clean[mark_name])
        results.append({
            "mark": mark_name,
            "pearson_r": round(r_p, 4),
            "pearson_p": f"{p_p:.2e}",
            "spearman_r": round(r_s, 4),
            "spearman_p": f"{p_s:.2e}",
        })
        print(f"  {mark_name}: Pearson r={r_p:.3f}, Spearman r={r_s:.3f}")

    corr_df = pd.DataFrame(results)
    corr_path = TABLES_DIR / f"{SAMPLE}_methylation_chipseq_correlations.tsv"
    corr_df.to_csv(corr_path, sep="\t", index=False)
    print(f"  Корреляции сохранены: {corr_path}")

    # Scatter plots
    fig, axes = plt.subplots(1, len(CHIP_TRACKS), figsize=(12, 5))
    if len(CHIP_TRACKS) == 1:
        axes = [axes]

    df_plot = df_clean.sample(min(50_000, len(df_clean)), random_state=42)

    for ax, (mark_name, _) in zip(axes, CHIP_TRACKS.items()):
        ax.scatter(
            df_plot["beta_value"], df_plot[mark_name],
            alpha=0.05, s=1, color="steelblue"
        )
        r = corr_df[corr_df["mark"] == mark_name]["pearson_r"].values[0]
        ax.set_xlabel("DNA methylation (beta-value)")
        ax.set_ylabel(f"{mark_name} signal (FE)")
        ax.set_title(f"{mark_name} vs methylation\nPearson r = {r}")

    plt.tight_layout()
    fig_path = FIGURES_DIR / f"{SAMPLE}_methylation_vs_chipseq.png"
    plt.savefig(fig_path, dpi=150)
    plt.close()
    print(f"  График сохранён: {fig_path}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("День 4: WGBS — анализ метилирования")
    print("=" * 60)

    filtered = run_part1()
    run_part2(filtered)
    run_part3()

    print("\n" + "=" * 60)
    print("✓ Пайплайн завершён!")
    print("=" * 60)
    print("\nФайлы для IGV:")
    for f in sorted(TRACKS_DIR.glob("*.bw")):
        print(f"  {f}")
    print("\nТаблицы:")
    for f in sorted((TABLES_DIR).glob("*")):
        print(f"  {f}")
    print("\nГрафики:")
    for f in sorted(FIGURES_DIR.glob("*.png")):
        print(f"  {f}")
