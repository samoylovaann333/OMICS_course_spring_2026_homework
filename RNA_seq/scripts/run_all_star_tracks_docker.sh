#!/usr/bin/env bash
set -euo pipefail

# STAR через Docker (native macOS STAR не читает FASTQ на этой машине — см. README).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THREADS="${THREADS:-4}"

GENOME_FA="${GENOME_FA:-${ROOT_DIR}/../day1_HiC_practice/data/reference/T2T_human.fna}"
CHROM_SIZES="${CHROM_SIZES:-${ROOT_DIR}/../day1_HiC_practice/data/reference/chrom.sizes}"
STAR_INDEX="${STAR_INDEX:-${ROOT_DIR}/../day1_HiC_practice/data/reference/star_index_sparse4}"
STAR_IMAGE="quay.io/biocontainers/star:2.7.11b--h43eeafb_2"

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

check_command docker
check_command samtools
check_command bedtools
check_command bedGraphToBigWig

run_star_docker() {
  local out_dir="$1"; shift
  local raw_dir="$1"; shift
  local prefix_name="$1"; shift
  docker run --rm \
    -v "${STAR_INDEX}:/genome:ro" \
    -v "${raw_dir}:/data:ro" \
    -v "${out_dir}:/out" \
    --entrypoint /usr/local/bin/STAR-plain \
    "${STAR_IMAGE}" \
    --runThreadN "${THREADS}" \
    --genomeDir /genome \
    --readFilesCommand zcat \
    --outFileNamePrefix "/out/${prefix_name}_" \
    "$@"
}

make_tracks() {
  local bam="$1"
  local assay="$2"
  local sample="$3"

  local bedgraph="${TRACK_DIR}/${assay}/${sample}.${assay}.STAR.bedGraph"
  local sorted_bedgraph="${TRACK_DIR}/${assay}/${sample}.${assay}.STAR.sorted.bedGraph"
  local bigwig="${TRACK_DIR}/${assay}/${sample}.${assay}.STAR.bw"

  bedtools genomecov -ibam "${bam}" -bg -split > "${bedgraph}"
  sort -k1,1 -k2,2n "${bedgraph}" > "${sorted_bedgraph}"
  bedGraphToBigWig "${sorted_bedgraph}" "${CHROM_SIZES}" "${bigwig}"
}

align_rnaseq() {
  local sample="$1"
  local r1="${RNA_RAW_DIR}/${sample}_R1.fastq.gz"
  local r2="${RNA_RAW_DIR}/${sample}_R2.fastq.gz"
  local out_dir="${STAR_DIR}/rnaseq"
  local bam="${out_dir}/${sample}.rnaseq.STAR.bam"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "Skipping RNA-seq ${sample}: FASTQ files not found."
    return
  fi

  if [[ -s "${bam}" && "${FORCE:-0}" != "1" ]]; then
    echo "RNA-seq ${sample}: BAM already exists. Skipping alignment."
  else
    run_star_docker "${out_dir}" "${RNA_RAW_DIR}" "${sample}" \
      --readFilesIn "/data/$(basename "$r1")" "/data/$(basename "$r2")" \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattrRGline "ID:${sample}_rnaseq" "SM:${sample}" "PL:ILLUMINA"

    mv -f "${out_dir}/${sample}_Aligned.sortedByCoord.out.bam" "${bam}"
    samtools index -@ "${THREADS}" "${bam}"
  fi

  samtools view -@ "${THREADS}" -b -q 30 "${bam}" \
    > "${out_dir}/${sample}.rnaseq.STAR.q30.bam"
  samtools index -@ "${THREADS}" "${out_dir}/${sample}.rnaseq.STAR.q30.bam"

  make_tracks "${bam}" "rnaseq" "${sample}"
}

align_cage() {
  local sample="$1"
  local r1="${CAGE_RAW_DIR}/${sample}_R1.fastq.gz"
  local out_dir="${STAR_DIR}/cage"
  local bam="${out_dir}/${sample}.cage.STAR.bam"

  if [[ ! -s "${r1}" ]]; then
    echo "Skipping CAGE ${sample}: FASTQ file not found."
    return
  fi

  if [[ -s "${bam}" && "${FORCE:-0}" != "1" ]]; then
    echo "CAGE ${sample}: BAM already exists. Skipping alignment."
  else
    run_star_docker "${out_dir}" "${CAGE_RAW_DIR}" "${sample}" \
      --readFilesIn "/data/$(basename "$r1")" \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattrRGline "ID:${sample}_cage" "SM:${sample}" "PL:ILLUMINA"

    mv -f "${out_dir}/${sample}_Aligned.sortedByCoord.out.bam" "${bam}"
    samtools index -@ "${THREADS}" "${bam}"
  fi

  samtools view -@ "${THREADS}" -b -q 30 "${bam}" \
    > "${out_dir}/${sample}.cage.STAR.q30.bam"
  samtools index -@ "${THREADS}" "${out_dir}/${sample}.cage.STAR.q30.bam"

  make_tracks "${bam}" "cage" "${sample}"
}

for sample in "${SAMPLES[@]}"; do
  echo "== ${sample}: RNA-seq =="
  align_rnaseq "${sample}"

  echo "== ${sample}: CAGE =="
  align_cage "${sample}"
done

echo "Done."
