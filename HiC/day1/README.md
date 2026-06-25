# Day 1 — Hi-C: подготовка ридов и получение `.hic` файлов

Задание: https://github.com/dpanc2/OMICS_course_spring_2026/tree/main/day1_HiC_practice

Референс — T2T-CHM13v2.0 (NCBI), хромосомы переименованы в `chr*`, плюс bwa-индекс,
`chrom.sizes` и сайты рестрикции DpnII. Дальше клонировала форк Juicer
(`juicer_course_version`) и прогнала 4 образца (`MoPh7`, `MoPh11`, `MoPh14`, `MoPh15`)
через [`scripts/run_hic_pipeline.sh`](scripts/run_hic_pipeline.sh): FastQC → cutadapt →
Juicer. Готовые `.hic` лежат в [`results/hic/`](results/hic).

## macOS-грабли

Juicer написан под gawk и GNU-утилиты, а на macOS по умолчанию BSD-версии:

1. Системный `awk` не понимает `and()`/`or()` из gawk — ставила `gawk` через brew
   и подкладывала его в `PATH` раньше системного (см. начало `run_hic_pipeline.sh`).
2. Системный `zcat` падает на обычных `.gz` без расширения `.Z`
   (`can't stat: file.gz (file.gz.Z): No such file or directory`) — добавила шим
   `tools/bin/zcat`, который просто вызывает `gzip -dc`.

## Как повторить

Сырые данные (геном, FASTQ, промежуточные файлы Juicer) в репозиторий не входят —
слишком тяжёлые.

```bash
conda create -n hic_practice -c conda-forge -c bioconda \
  fastqc cutadapt bwa samtools openjdk=11 wget
conda activate hic_practice
brew install gawk

mkdir -p data/raw data/trimmed data/reference data/juicer
mkdir -p results/fastqc_raw results/cutadapt results/hic

# геном T2T-CHM13v2.0
wget -O data/reference/T2T_human.fna.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
gzip -dkf data/reference/T2T_human.fna.gz
python3 scripts/rename_chroms_t2t.py
bwa index data/reference/T2T_human.fna
samtools faidx data/reference/T2T_human.fna
cut -f1,2 data/reference/T2T_human.fna.fai > data/reference/chrom.sizes

# Juicer
git clone --branch juicer_course_version --single-branch \
  https://github.com/dpanc2/OMICS_course_spring_2026.git tools/juicer
python3 tools/juicer/misc/generate_site_positions.py DpnII T2T_human data/reference/T2T_human.fna
mv T2T_human_DpnII.txt data/reference/restriction_sites_DpnII.txt

# zcat-шим, только для macOS
mkdir -p tools/bin
printf '#!/usr/bin/env bash\nexec gzip -dc "$@"\n' > tools/bin/zcat
chmod +x tools/bin/zcat

bash scripts/run_hic_pipeline.sh MoPh7 MoPh11 MoPh14 MoPh15
```

## По итоговому вопросу задания

Смотрела на 4 карты в Juicebox — но в учебных FASTQ всего ~45–65 тыс. ридов на образец
(вместо обычных десятков миллионов), поэтому карты получаются очень разреженными и
по ним нельзя надёжно судить о крупных структурных перестройках. Для этого нужны
полноразмерные данные — как в [`../day2`](../day2), где используются уже готовые
"enriched" карты с 10–14 млн контактов на образец.
