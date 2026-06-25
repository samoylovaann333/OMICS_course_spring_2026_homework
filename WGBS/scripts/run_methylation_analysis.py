#!/usr/bin/env python3
"""
День 4. Анализ метилирования
Запуск: python3 run_methylation_analysis.py
"""

import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# ПУТИ
# =============================================================================

ROOT = Path(__file__).parent
REF_DIR = ROOT.parent / "day1_HiC_practice" / "data" / "reference"

GTF_FILE = REF_DIR / "chm13v2.0.gff3.gz"
METH_TABLE = ROOT / "results" / "tables" / "MoPh7_cpg_methylation_values.tsv.gz"

FIGURES_DIR = ROOT / "results" / "figures"
TABLES_DIR = ROOT / "results" / "tables"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE = "MoPh7"
PROMOTER_UP = 1000
MIN_CPG = 3

# =============================================================================
# ФУНКЦИИ
# =============================================================================

def get_promoters_and_genes(gff3_path):
    """Извлекает промоторы и тела генов из GFF3 файла."""
    print("📂 Читаем GFF3 аннотацию...")
    
    promoters = []
    gene_bodies = []
    
    opener = gzip.open if str(gff3_path).endswith(".gz") else open
    
    with opener(gff3_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            
            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            
            # Извлекаем имя гена
            gene_name = ""
            for attr in attrs.split(";"):
                if attr.startswith("Name="):
                    gene_name = attr.split("=")[1]
                    break
            
            # Промотор (1000 bp перед TSS)
            if strand == "+":
                tss = start
                prom_start = max(0, tss - PROMOTER_UP)
                prom_end = tss
                body_start = start + 1000
                body_end = end
            else:
                tss = end
                prom_start = max(0, tss - PROMOTER_UP)
                prom_end = tss
                body_start = start
                body_end = end - 1000
            
            # Добавляем промотор
            if prom_start < prom_end:
                promoters.append({
                    "chrom": chrom,
                    "start": prom_start,
                    "end": prom_end,
                    "gene_name": gene_name,
                    "strand": strand,
                    "tss": tss
                })
            
            # Добавляем тело гена
            if body_start < body_end:
                gene_bodies.append({
                    "chrom": chrom,
                    "start": body_start,
                    "end": body_end
                })
    
    promoters_df = pd.DataFrame(promoters)
    gene_bodies_df = pd.DataFrame(gene_bodies)
    
    print(f"  ✅ Найдено генов: {len(promoters_df):,}")
    print(f"  ✅ Тел генов: {len(gene_bodies_df):,}")
    
    return promoters_df, gene_bodies_df


def load_methylation_data(path):
    """Загружает таблицу метилирования."""
    print("📂 Читаем данные метилирования...")
    df = pd.read_csv(path, sep="\t")
    print(f"  ✅ CpG сайтов: {len(df):,}")
    return df


def calculate_region_methylation(meth_df, regions_df, min_cpg=3):
    """Считает среднее метилирование для каждого региона."""
    print("  ⏳ Считаем метилирование в регионах...")
    
    # Группируем CpG по хромосомам
    meth_by_chrom = {chrom: grp for chrom, grp in meth_df.groupby("chrom")}
    
    results = []
    
    for idx, region in regions_df.iterrows():
        chrom = region["chrom"]
        
        if chrom not in meth_by_chrom:
            continue
        
        cpg = meth_by_chrom[chrom]
        mask = (cpg["start"] >= region["start"]) & (cpg["end"] <= region["end"])
        cpg_in_region = cpg[mask]
        
        if len(cpg_in_region) < min_cpg:
            continue
        
        row = region.to_dict()
        row["mean_beta"] = cpg_in_region["beta_value"].mean()
        row["n_cpg"] = len(cpg_in_region)
        results.append(row)
    
    result_df = pd.DataFrame(results)
    print(f"  ✅ Регионов с данными: {len(result_df):,}")
    return result_df


def get_background_methylation(meth_df, promoters_df, gene_bodies_df, n_sample=100000):
    """Выбирает случайные CpG вне генов (фон)."""
    print("  ⏳ Выбираем фоновые CpG...")
    
    # Создаем множество занятых позиций
    occupied = set()
    
    for _, row in promoters_df.iterrows():
        for pos in range(row["start"], row["end"] + 1, 100):
            occupied.add((row["chrom"], pos))
    
    for _, row in gene_bodies_df.iterrows():
        for pos in range(row["start"], row["end"] + 1, 100):
            occupied.add((row["chrom"], pos))
    
    # Фильтруем CpG, не попадающие в занятые позиции
    def is_background(row):
        return (row["chrom"], row["start"]) not in occupied
    
    bg_mask = meth_df.apply(is_background, axis=1)
    bg = meth_df[bg_mask]
    
    if len(bg) > n_sample:
        bg = bg.sample(n_sample, random_state=42)
    elif len(bg) == 0:
        print("  ⚠️ Не найдено фоновых CpG! Использую случайные из всех.")
        bg = meth_df.sample(min(n_sample, len(meth_df)), random_state=42)
    
    print(f"  ✅ Фоновых CpG: {len(bg):,}")
    return bg["beta_value"]


def plot_methylation_comparison(prom_beta, body_beta, bg_beta, out_dir):
    """Строит графики сравнения метилирования."""
    print("  ⏳ Строим графики...")
    
    # Собираем данные
    data = pd.DataFrame({
        "beta_value": list(prom_beta) + list(body_beta) + list(bg_beta),
        "region": ["Промоторы"] * len(prom_beta) + \
                  ["Тела генов"] * len(body_beta) + \
                  ["Фон"] * len(bg_beta)
    })
    
    data = data.dropna()
    
    if len(data) == 0:
        print("  ❌ Нет данных для построения графиков!")
        return
    
    order = ["Промоторы", "Тела генов", "Фон"]
    colors = ["#4C9BE8", "#E8844C", "#8FBE5A"]
    
    # График 1: Полный диапазон (0-1)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=data, x="region", y="beta_value",
                order=order, palette=colors, ax=ax, showfliers=False)
    ax.set_title(f"{SAMPLE}: Метилирование по типам регионов", fontsize=14)
    ax.set_ylabel("Beta-value (метилирование)", fontsize=12)
    ax.set_ylim(0, 1)
    
    # Добавляем медианы
    for i, region in enumerate(order):
        region_data = data[data["region"] == region]["beta_value"]
        if len(region_data) > 0:
            med = region_data.median()
            ax.text(i, med + 0.02, f"{med:.2f}", ha="center", 
                   fontsize=11, fontweight="bold")
    
    plt.tight_layout()
    plt.savefig(out_dir / f"{SAMPLE}_methylation_full.png", dpi=150)
    plt.close()
    
    # График 2: Zoom (0.5-1.0)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=data, x="region", y="beta_value",
                order=order, palette=colors, ax=ax, showfliers=False)
    ax.set_title(f"{SAMPLE}: Метилирование (zoom 0.5-1.0)", fontsize=14)
    ax.set_ylabel("Beta-value (метилирование)", fontsize=12)
    ax.set_ylim(0.5, 1.0)
    ax.axhline(y=0.8, color='gray', linestyle='--', alpha=0.5)
    
    for i, region in enumerate(order):
        region_data = data[data["region"] == region]["beta_value"]
        if len(region_data) > 0:
            med = region_data.median()
            ax.text(i, med + 0.005, f"{med:.2f}", ha="center",
                   fontsize=11, fontweight="bold")
    
    plt.tight_layout()
    plt.savefig(out_dir / f"{SAMPLE}_methylation_zoom.png", dpi=150)
    plt.close()
    
    # График 3: Violin plot
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.violinplot(data=data, x="region", y="beta_value",
                   order=order, palette=colors, ax=ax, inner="quartile", cut=0)
    ax.set_title(f"{SAMPLE}: Распределение метилирования", fontsize=14)
    ax.set_ylabel("Beta-value (метилирование)", fontsize=12)
    ax.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig(out_dir / f"{SAMPLE}_methylation_violin.png", dpi=150)
    plt.close()
    
    print(f"  ✅ Графики сохранены в {out_dir}")


def plot_tss_profile(meth_df, promoters_df, out_dir):
    """Строит профиль метилирования вокруг TSS."""
    print("  ⏳ Строим TSS профиль...")
    
    # Берем до 5000 генов
    n_genes = min(5000, len(promoters_df))
    genes = promoters_df.sample(n_genes, random_state=42)
    
    window = 2000
    n_bins = 40
    bin_size = (2 * window) // n_bins
    
    # Группируем CpG по хромосомам
    meth_by_chrom = {chrom: grp.set_index("start") 
                     for chrom, grp in meth_df.groupby("chrom")}
    
    profiles = []
    
    for _, gene in genes.iterrows():
        chrom = gene["chrom"]
        tss = gene["tss"]
        strand = gene["strand"]
        
        if chrom not in meth_by_chrom:
            continue
        
        cpg = meth_by_chrom[chrom]
        
        # Берем CpG в окне вокруг TSS
        region_cpg = cpg[
            (cpg.index >= tss - window) & (cpg.index <= tss + window)
        ].copy()
        
        if len(region_cpg) < 5:
            continue
        
        region_cpg = region_cpg.reset_index()
        region_cpg["rel_pos"] = region_cpg["start"] - tss
        
        # Для генов на минус-цепи переворачиваем
        if strand == "-":
            region_cpg["rel_pos"] = -region_cpg["rel_pos"]
        
        # Биннируем
        bins = np.linspace(-window, window, n_bins + 1)
        region_cpg["bin"] = pd.cut(region_cpg["rel_pos"], bins=bins, labels=False)
        bin_means = region_cpg.groupby("bin")["beta_value"].mean()
        
        profile = np.full(n_bins, np.nan)
        for b, val in bin_means.items():
            if pd.notna(b):
                profile[int(b)] = val
        
        profiles.append(profile)
    
    if not profiles:
        print("  ❌ Недостаточно данных для TSS профиля")
        return
    
    profiles = np.array(profiles)
    mean_profile = np.nanmean(profiles, axis=0)
    sem_profile = np.nanstd(profiles, axis=0) / np.sqrt(
        np.sum(~np.isnan(profiles), axis=0)
    )
    
    bin_centers = np.linspace(-window + bin_size/2, window - bin_size/2, n_bins)
    
    # Строим график
    fig, ax = plt.subplots(figsize=(10, 5))
    
    ax.plot(bin_centers, mean_profile, color="#4C9BE8", 
            linewidth=2.5, label="Среднее метилирование")
    ax.fill_between(
        bin_centers,
        mean_profile - sem_profile,
        mean_profile + sem_profile,
        alpha=0.25, color="#4C9BE8", label="±SEM"
    )
    
    ax.axvline(0, color="red", linestyle="--", linewidth=1.5, label="TSS")
    ax.set_xlabel("Позиция относительно TSS (bp)", fontsize=12)
    ax.set_ylabel("Среднее метилирование (beta-value)", fontsize=12)
    ax.set_title(f"{SAMPLE}: Профиль метилирования вокруг TSS", fontsize=13)
    ax.legend()
    ax.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig(out_dir / f"{SAMPLE}_tss_profile.png", dpi=150)
    plt.close()
    
    print(f"  ✅ TSS профиль сохранен")


# =============================================================================
# ОСНОВНАЯ ПРОГРАММА
# =============================================================================

def main():
    print("=" * 60)
    print("ДЕНЬ 4: АНАЛИЗ МЕТИЛИРОВАНИЯ")
    print("=" * 60)
    
    # 1. Загружаем данные
    print("\n[1/5] Загрузка данных...")
    meth_df = load_methylation_data(METH_TABLE)
    
    # 2. Извлекаем промоторы и тела генов
    print("\n[2/5] Извлечение промоторов и тел генов...")
    promoters_df, gene_bodies_df = get_promoters_and_genes(GTF_FILE)
    
    # Оставляем только основные хромосомы
    valid_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    promoters_df = promoters_df[promoters_df["chrom"].isin(valid_chroms)]
    gene_bodies_df = gene_bodies_df[gene_bodies_df["chrom"].isin(valid_chroms)]
    
    print(f"  ✅ После фильтрации: {len(promoters_df):,} промоторов")
    
    # 3. Считаем метилирование промоторов
    print("\n[3/5] Расчет метилирования промоторов...")
    promoter_meth = calculate_region_methylation(meth_df, promoters_df)
    
    if len(promoter_meth) > 0:
        median_prom = promoter_meth["mean_beta"].median()
        print(f"  📊 Медиана метилирования промоторов: {median_prom:.3f}")
        
        # Топ-10 неметилированных промоторов
        top_low = promoter_meth.nsmallest(10, "mean_beta")
        print("\n  🏆 Топ-10 наименее метилированных промоторов:")
        print(top_low[["gene_name", "chrom", "start", "end", "mean_beta", "n_cpg"]].to_string(index=False))
        
        # Сохраняем таблицу
        out_table = TABLES_DIR / f"{SAMPLE}_promoter_methylation.tsv"
        promoter_meth.to_csv(out_table, sep="\t", index=False)
        print(f"\n  ✅ Таблица сохранена: {out_table}")
    
    # 4. Сравниваем с телами генов и фоном
    print("\n[4/5] Сравнение с телами генов и фоном...")
    
    # Тела генов
    gene_meth = calculate_region_methylation(meth_df, gene_bodies_df)
    median_gene = gene_meth["mean_beta"].median() if len(gene_meth) > 0 else np.nan
    print(f"  📊 Медиана тел генов: {median_gene:.3f} (n={len(gene_meth):,})")
    
    # Фон
    bg_beta = get_background_methylation(meth_df, promoters_df, gene_bodies_df)
    median_bg = bg_beta.median()
    print(f"  📊 Медиана фона: {median_bg:.3f} (n={len(bg_beta):,})")
    
    # Статистический тест
    if len(promoter_meth) > 0 and len(gene_meth) > 0 and len(bg_beta) > 0:
        try:
            stat, pval = stats.kruskal(
                promoter_meth["mean_beta"],
                gene_meth["mean_beta"],
                bg_beta
            )
            print(f"\n  📊 Kruskal-Wallis тест:")
            print(f"     H = {stat:.1f}, p = {pval:.2e}")
            if pval < 0.001:
                print("     ✅ Различия статистически значимы (p < 0.001)")
        except Exception as e:
            print(f"  ⚠️ Ошибка в статистическом тесте: {e}")
    
    # 5. Строим графики
    print("\n[5/5] Построение графиков...")
    
    if len(promoter_meth) > 0 and len(gene_meth) > 0 and len(bg_beta) > 0:
        # Графики сравнения
        plot_methylation_comparison(
            promoter_meth["mean_beta"],
            gene_meth["mean_beta"],
            bg_beta,
            FIGURES_DIR
        )
        
        # TSS профиль
        plot_tss_profile(meth_df, promoters_df, FIGURES_DIR)
    
    # Итог
    print("\n" + "=" * 60)
    print("✅ АНАЛИЗ ЗАВЕРШЕН!")
    print("=" * 60)
    print("\n📁 Результаты:")
    print(f"   Графики: {FIGURES_DIR}")
    print(f"   Таблицы: {TABLES_DIR}")
    print("\n📊 Файлы:")
    for f in sorted(FIGURES_DIR.glob("*.png")):
        print(f"     - {f.name}")
    for f in sorted(TABLES_DIR.glob("*.tsv")):
        print(f"     - {f.name}")

if __name__ == "__main__":
    main()
