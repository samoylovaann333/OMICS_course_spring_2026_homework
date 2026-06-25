# День 4: WGBS — анализ метилирования ДНК

## Структура

```
day4/
├── notebooks/
│   ├── bismark_qc_and_methylation_values.ipynb
│   ├── create_bigwig_tracks_for_igv.ipynb
│   └── integrate_methylation_with_h3k27ac.ipynb
├── results/
│   ├── MoPh7_methylation_distributions.png
│   ├── MoPh7_methylation_full.png
│   ├── MoPh7_methylation_violin.png
│   ├── MoPh7_methylation_zoom.png
│   ├── MoPh7_promoter_methylation.tsv
│   └── MoPh7_tss_profile.png
└── scripts/
    ├── run_pipeline.py
    └── run_methylation_analysis.py
```

## Данные

Входные данные — готовый output Bismark для образца MoPh7. Скачиваются с FTP-сервера:

```bash
python3 scripts/download_bismark_data.py
```

Референсный геном: T2T-CHM13v2.0 (`chm13v2.0.fa`), аннотация: `chm13v2.0.gff3.gz`

## Скрипты

`run_pipeline.py` — основной пайплайн, объединяет логику всех трёх ноутбуков:
- Читает `bismark.cov.gz`, считает beta-value и M-value для каждого CpG
- Фильтрует CpG с coverage < 5
- Строит QC-графики распределений
- Создаёт bigWig треки для IGV (beta, M-value, coverage, GC-content, CpG obs/exp)

`run_methylation_analysis.py` — дополнительный анализ:
- Извлекает промоторы из GFF3 (1000 bp перед TSS с учётом ориентации гена)
- Считает среднее метилирование для 55 981 промоторов (медиана = 0.61)
- Сравнивает метилирование промоторов, тел генов и фона
- Строит метапрофиль метилирования вокруг TSS (±2000 bp, 5000 генов)

## Запуск

```bash
conda activate hic_practice
python3 scripts/run_pipeline.py
python3 scripts/run_methylation_analysis.py
```

## Окружение

```bash
conda create -n hic_practice python=3.9 -c conda-forge -y
conda activate hic_practice
conda install -c conda-forge -c bioconda \
  pandas numpy matplotlib seaborn scipy jupyterlab pybigwig pysam -y
```
