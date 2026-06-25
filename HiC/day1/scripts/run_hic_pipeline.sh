#!/usr/bin/env bash
set -euo pipefail

# Hi-C pipeline: сырые paired-end риды -> .hic файл (Juicer)
# Использование: ./scripts/run_hic_pipeline.sh sample1 sample2 ...
# Для каждого <sample> ожидается, что в RAW_BASE_URL лежат файлы:
#   "Copy of <sample_tag>_R1_001.fastq.gz" / "Copy of <sample_tag>_R2_001.fastq.gz"
# Соответствие sample -> sample_tag берется из SAMPLE_TAGS ниже.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

# Juicer-скрипты используют gawk-расширения (and(), or()), которых нет в
# системном BSD awk на macOS, поэтому подкладываем gawk впереди в PATH.
# Также подменяем системный /usr/bin/zcat (который на macOS не умеет читать
# .gz без расширения .Z) на shim вокруг "gzip -dc".
export PATH="$ROOT_DIR/tools/bin:/opt/homebrew/opt/gawk/libexec/gnubin:$PATH"

RAW_BASE_URL="https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics"

GENOME_NAME="T2T_human"
GENOME_FASTA="$ROOT_DIR/data/reference/T2T_human.fna"
CHROM_SIZES="$ROOT_DIR/data/reference/chrom.sizes"
RESTRICTION_SITES="$ROOT_DIR/data/reference/restriction_sites_DpnII.txt"
ENZYME="DpnII"
JUICER_DIR="$ROOT_DIR/tools/juicer"
THREADS=4
CUTADAPT_BIN="/opt/homebrew/Caskroom/miniforge/base/envs/hic_practice/bin/cutadapt"

# sample -> исходный тег файла на сервере (часть имени до _R1/_R2_001.fastq.gz)
sample_tag() {
  case "$1" in
    MoPh7)  echo "MoPh7_S85_L001" ;;
    MoPh11) echo "MoPh11_S86_L001" ;;
    MoPh14) echo "MoPh14_S87_L001" ;;
    MoPh15) echo "MoPh15_S88_L001" ;;
    *) echo "ERROR: неизвестный sample '$1', добавьте его тег в sample_tag()" >&2; exit 1 ;;
  esac
}

if [ "$#" -eq 0 ]; then
  echo "Использование: $0 <sample1> [sample2 ...]" >&2
  exit 1
fi

for sample in "$@"; do
  tag="$(sample_tag "$sample")"
  echo "=== [$sample] Шаг 1/6: скачивание FASTQ ==="
  mkdir -p data/raw
  raw_r1="data/raw/${sample}_R1.fastq.gz"
  raw_r2="data/raw/${sample}_R2.fastq.gz"
  if [ ! -s "$raw_r1" ]; then
    curl -k -s -o "$raw_r1" "${RAW_BASE_URL}/Copy%20of%20${tag}_R1_001.fastq.gz"
  fi
  if [ ! -s "$raw_r2" ]; then
    curl -k -s -o "$raw_r2" "${RAW_BASE_URL}/Copy%20of%20${tag}_R2_001.fastq.gz"
  fi

  echo "=== [$sample] Шаг 2/6: FastQC на сырых ридах ==="
  mkdir -p results/fastqc_raw
  fastqc "$raw_r1" "$raw_r2" -o results/fastqc_raw

  echo "=== [$sample] Шаг 3/6: обрезка адаптеров cutadapt ==="
  mkdir -p data/trimmed results/cutadapt
  trimmed_r1="data/trimmed/${sample}_R1.trimmed.fastq.gz"
  trimmed_r2="data/trimmed/${sample}_R2.trimmed.fastq.gz"
  "$CUTADAPT_BIN" \
    -q 20 \
    -m 70 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -o "$trimmed_r1" \
    -p "$trimmed_r2" \
    "$raw_r1" "$raw_r2" \
    > "results/cutadapt/${sample}.cutadapt.log" 2>&1

  echo "=== [$sample] Шаг 4/6: подготовка директории Juicer ==="
  juicer_fastq_dir="data/juicer/${sample}/fastq"
  mkdir -p "$juicer_fastq_dir"
  ln -sf "$ROOT_DIR/$trimmed_r1" "$juicer_fastq_dir/${sample}_R1.fastq.gz"
  ln -sf "$ROOT_DIR/$trimmed_r2" "$juicer_fastq_dir/${sample}_R2.fastq.gz"

  echo "=== [$sample] Шаг 5/6: запуск Juicer ==="
  bash "$JUICER_DIR/scripts/juicer.sh" \
    -D "$JUICER_DIR" \
    -d "$ROOT_DIR/data/juicer/${sample}" \
    -g "$GENOME_NAME" \
    -z "$GENOME_FASTA" \
    -p "$CHROM_SIZES" \
    -y "$RESTRICTION_SITES" \
    -s "$ENZYME" \
    -t "$THREADS"

  echo "=== [$sample] Шаг 6/6: сохранение .hic в results/hic ==="
  mkdir -p results/hic
  cp "data/juicer/${sample}/aligned/inter_30.hic" "results/hic/${sample}.inter_30.hic"

  echo "=== [$sample] готово: results/hic/${sample}.inter_30.hic ==="
done

echo "Все образцы обработаны: $*"
