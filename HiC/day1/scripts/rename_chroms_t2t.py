#!/usr/bin/env python3
"""Rename NCBI T2T-CHM13 FASTA records to chr-style chromosome names."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


DEFAULT_FASTA = Path("data/reference/T2T_human.fna")
DEFAULT_MAP = Path("data/reference/T2T_human.rename_chroms.tsv")


def chr_name_from_header(header: str) -> str | None:
    """Return chr-style name inferred from an NCBI FASTA header."""
    record_id, _, description = header.partition(" ")

    if record_id.startswith("chr"):
        return record_id

    chromosome_match = re.search(
        r"\bchromosome\s+([0-9]+|X|Y|M|MT)\b",
        description,
        flags=re.IGNORECASE,
    )
    if chromosome_match:
        chromosome = chromosome_match.group(1).upper()
        if chromosome == "MT":
            chromosome = "M"
        return f"chr{chromosome}"

    if re.search(r"\bmitochondr", description, flags=re.IGNORECASE):
        return "chrM"

    return None


def rename_fasta_in_place(fasta_path: Path, map_path: Path) -> None:
    tmp_path = fasta_path.with_suffix(fasta_path.suffix + ".tmp")
    seen_names: set[str] = set()

    with fasta_path.open() as source, tmp_path.open("w") as target, map_path.open("w") as mapping:
        mapping.write("old_id\tnew_id\told_header\n")

        for line in source:
            if not line.startswith(">"):
                target.write(line)
                continue

            old_header = line[1:].rstrip("\n")
            old_id = old_header.split()[0]
            new_id = chr_name_from_header(old_header)

            if new_id is None:
                raise ValueError(
                    "Cannot infer chr-style name for FASTA record: "
                    f"{old_header!r}"
                )

            if new_id in seen_names and new_id != old_id:
                raise ValueError(f"Duplicate chromosome name after renaming: {new_id}")

            seen_names.add(new_id)
            target.write(f">{new_id}\n")
            mapping.write(f"{old_id}\t{new_id}\t{old_header}\n")

    tmp_path.replace(fasta_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rename T2T-CHM13 FASTA headers from NCBI IDs to chr-style names."
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=DEFAULT_FASTA,
        help=f"Input FASTA to rename in place. Default: {DEFAULT_FASTA}",
    )
    parser.add_argument(
        "--map",
        type=Path,
        default=DEFAULT_MAP,
        help=f"Output TSV with old and new names. Default: {DEFAULT_MAP}",
    )
    args = parser.parse_args()

    if not args.fasta.exists():
        raise FileNotFoundError(f"FASTA file not found: {args.fasta}")

    args.map.parent.mkdir(parents=True, exist_ok=True)
    rename_fasta_in_place(args.fasta, args.map)
    print(f"Renamed FASTA headers in: {args.fasta}")
    print(f"Saved chromosome-name map to: {args.map}")


if __name__ == "__main__":
    main()
