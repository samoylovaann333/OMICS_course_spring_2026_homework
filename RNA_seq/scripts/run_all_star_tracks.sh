#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THREADS="${THREADS:-8}"

GENOME_FA="${GENOME_FA:-${ROOT_DIR}/../day1_HiC_practice/data/reference/T2T_human.fna}"
CHROM_SIZES="${CHROM_SIZES:-${ROOT_DIR}/../day1_HiC_practice/data/reference/chrom.sizes}"
STAR_INDEX="${STAR_INDEX:-${ROOT_DIR}/../day1_HiC_practice/data/reference/star_index}"

RNA_RAW_DIR="${ROOT_DIR}/data/raw/rnaseq"
CAGE_RAW_DIR="${ROOT_DIR}/data/raw/cage"

STAR_DIR="${ROOT_DIR}/results/star"
TRACK_DIR="${ROOT_DIR}/results/tracks"
LOG_DIR="${ROOT_DIR}/results/logs"

SAMPLES=(MoPh7 MoPh11 MoPh14 MoPh15)

mkdir -p "${STAR_DIR}/rnaseq" "${STAR_DIR}/cage" "${TRACK_DIR}/rnaseq" "${TRACK_DIR}/cage" "${LOG_DIR}"

check_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: command not found: $1" >&2
    exit 1
  fi
}

check_command STAR
check_command samtools
check_command bedtools
check_command bedGraphToBigWig

if [[ ! -d "${STAR_INDEX}" ]]; then
  echo "ERROR: STAR index not found: ${STAR_INDEX}" >&2
  echo "Build it first using the README instructions." >&2
  exit 1
fi

if [[ ! -s "${CHROM_SIZES}" ]]; then
  echo "ERROR: chromosome sizes file not found: ${CHROM_SIZES}" >&2
  exit 1
fi

make_tracks() {
  local bam="$1"
  local assay="$2"
  local sample="$3"

  local bedgraph="${TRACK_DIR}/${assay}/${sample}.${assay}.STAR.bedGraph"
  local sorted_bedgraph="${TRACK_DIR}/${assay}/${sample}.${assay}.STAR.sorted.bedGraph"
  local bigwig="${TRACK_DIR}/${assay}/${sample}.${assay}.STAR.bw"

  bedtools genomecov \
    -ibam "${bam}" \
    -bg \
    -split \
    > "${bedgraph}"

  sort -k1,1 -k2,2n "${bedgraph}" > "${sorted_bedgraph}"
  bedGraphToBigWig "${sorted_bedgraph}" "${CHROM_SIZES}" "${bigwig}"
}

align_rnaseq() {
  local sample="$1"
  local r1="${RNA_RAW_DIR}/${sample}_R1.fastq.gz"
  local r2="${RNA_RAW_DIR}/${sample}_R2.fastq.gz"
  local prefix="${STAR_DIR}/rnaseq/${sample}_"
  local bam="${STAR_DIR}/rnaseq/${sample}.rnaseq.STAR.bam"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "Skipping RNA-seq ${sample}: FASTQ files not found."
    return
  fi

  if [[ -s "${bam}" && "${FORCE:-0}" != "1" ]]; then
    echo "RNA-seq ${sample}: BAM already exists. Skipping alignment."
  else
    STAR \
      --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX}" \
      --readFilesIn "${r1}" "${r2}" \
      --readFilesCommand gzip -dc \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "${prefix}" \
      --outSAMattrRGline "ID:${sample}_rnaseq" "SM:${sample}" "PL:ILLUMINA"

    mv -f "${prefix}Aligned.sortedByCoord.out.bam" "${bam}"
    samtools index -@ "${THREADS}" "${bam}"
  fi

  samtools view -@ "${THREADS}" -b -q 30 "${bam}" \
    > "${STAR_DIR}/rnaseq/${sample}.rnaseq.STAR.q30.bam"
  samtools index -@ "${THREADS}" "${STAR_DIR}/rnaseq/${sample}.rnaseq.STAR.q30.bam"

  make_tracks "${bam}" "rnaseq" "${sample}"
}

align_cage() {
  local sample="$1"
  local r1="${CAGE_RAW_DIR}/${sample}_R1.fastq.gz"
  local prefix="${STAR_DIR}/cage/${sample}_"
  local bam="${STAR_DIR}/cage/${sample}.cage.STAR.bam"

  if [[ ! -s "${r1}" ]]; then
    echo "Skipping CAGE ${sample}: FASTQ file not found."
    return
  fi

  if [[ -s "${bam}" && "${FORCE:-0}" != "1" ]]; then
    echo "CAGE ${sample}: BAM already exists. Skipping alignment."
  else
    STAR \
      --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX}" \
      --readFilesIn "${r1}" \
      --readFilesCommand gzip -dc \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "${prefix}" \
      --outSAMattrRGline "ID:${sample}_cage" "SM:${sample}" "PL:ILLUMINA"

    mv -f "${prefix}Aligned.sortedByCoord.out.bam" "${bam}"
    samtools index -@ "${THREADS}" "${bam}"
  fi

  samtools view -@ "${THREADS}" -b -q 30 "${bam}" \
    > "${STAR_DIR}/cage/${sample}.cage.STAR.q30.bam"
  samtools index -@ "${THREADS}" "${STAR_DIR}/cage/${sample}.cage.STAR.q30.bam"

  make_tracks "${bam}" "cage" "${sample}"
}

for sample in "${SAMPLES[@]}"; do
  echo "== ${sample}: RNA-seq =="
  align_rnaseq "${sample}"

  echo "== ${sample}: CAGE =="
  align_cage "${sample}"
done

echo "Done."
