#!/usr/bin/env python3
"""Download RNA-seq and CAGE raw reads for the RNA practice."""

from __future__ import annotations

import argparse
from pathlib import Path

import requests
import urllib3
from tqdm import tqdm


VERIFY_SSL = False

RAW_FILES = {
    "rnaseq": {
        "MoPh7": {
            "R1": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2026-03-24-BGI-RNAseq/F25A430000584_HOMckpbR/7DR5/7DR5_1.fq.gz",
            "R2": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2026-03-24-BGI-RNAseq/F25A430000584_HOMckpbR/7DR5/7DR5_2.fq.gz",
        },
        "MoPh11": {
            "R1": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2026-03-24-BGI-RNAseq/F26A430000218_HOMynzpR/11TR6/11TR6_1.fq.gz",
            "R2": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2026-03-24-BGI-RNAseq/F26A430000218_HOMynzpR/11TR6/11TR6_2.fq.gz",
        },
        "MoPh14": {
            "R1": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2026-03-24-BGI-RNAseq/F25A430000584_HOMckpbR/14DR7/14DR7_1.fq.gz",
            "R2": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2026-03-24-BGI-RNAseq/F25A430000584_HOMckpbR/14DR7/14DR7_2.fq.gz",
        },
        "MoPh15": {
            "R1": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2026-03-24-BGI-RNAseq/F25A430000584_HOMckpbR/15DR8/15DR8_1.fq.gz",
            "R2": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2026-03-24-BGI-RNAseq/F25A430000584_HOMckpbR/15DR8/15DR8_2.fq.gz",
        },
    },
    "cage": {
        "MoPh7": {
            "R1": "https://genedev.bionet.nsc.ru/ftp/_RawReads/KAZAN_2025/Moph-7_S5_L001_R1_001.fastq.gz",
        },
        "MoPh11": {
            "R1": "https://genedev.bionet.nsc.ru/ftp/_RawReads/KAZAN_2025/Moph-11_S6_L001_R1_001.fastq.gz",
        },
        "MoPh14": {
            "R1": "https://genedev.bionet.nsc.ru/ftp/_RawReads/KAZAN_2025/Moph-14_S7_L001_R1_001.fastq.gz",
        },
        "MoPh15": {
            "R1": "https://genedev.bionet.nsc.ru/ftp/_RawReads/KAZAN_2025/Moph-15_S8_L001_R1_001.fastq.gz",
        },
    },
}


def download_file(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"{out_path} already exists. Skipping.")
        return

    print(f"Downloading {out_path.name}")
    response = requests.get(url, stream=True, verify=VERIFY_SSL)
    response.raise_for_status()

    total_size = int(response.headers.get("content-length", 0))
    block_size = 1024 * 1024

    with out_path.open("wb") as handle, tqdm(
        total=total_size,
        unit="iB",
        unit_scale=True,
        unit_divisor=1024,
        desc=out_path.name,
    ) as progress:
        for chunk in response.iter_content(block_size):
            if not chunk:
                continue
            handle.write(chunk)
            progress.update(len(chunk))


def selected_assays(value: str) -> list[str]:
    if value == "all":
        return list(RAW_FILES)
    return [value]


def selected_samples(value: str, assay: str) -> list[str]:
    if value == "all":
        return list(RAW_FILES[assay])
    if value not in RAW_FILES[assay]:
        available = ", ".join(RAW_FILES[assay])
        raise ValueError(f"Unknown sample {value!r} for {assay}. Available: {available}")
    return [value]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--assay",
        choices=["rnaseq", "cage", "all"],
        default="all",
        help="Which data type to download. Default: all.",
    )
    parser.add_argument(
        "--sample",
        default="all",
        help="Sample name: MoPh7, MoPh11, MoPh14, MoPh15, or all. Default: all.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data/raw"),
        help="Output directory. Default: data/raw.",
    )
    args = parser.parse_args()

    if not VERIFY_SSL:
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    for assay in selected_assays(args.assay):
        for sample in selected_samples(args.sample, assay):
            for read_name, url in RAW_FILES[assay][sample].items():
                out_name = f"{sample}_{read_name}.fastq.gz"
                out_path = args.data_dir / assay / out_name
                download_file(url, out_path)

    print("Done.")


if __name__ == "__main__":
    main()
