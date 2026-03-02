#!/usr/bin/env python3
"""
Produce a BED file of introns from an exon BED input.

Input format (tab-separated):
  chrom, start, end, gene, exon_num, strand, transcript_id, [extra_field...]

Introns are computed between consecutive exons within each transcript,
ordered by exon number.
"""

import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Produce intron BED from exon BED file"
    )
    parser.add_argument(
        "exon_bed",
        help="Input exon BED file (tab-separated)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output intron BED file (default: stdout)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Group exons by (chromosome, transcript_id, strand)
    groups = defaultdict(list)

    with open(args.exon_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 7:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            gene = fields[3]
            exon_num = int(fields[4])
            strand = fields[5]
            transcript_id = fields[6]
            key = (chrom, transcript_id, strand)
            groups[key].append((exon_num, start, end, gene))

    # Sort exons by exon number within each group and compute introns
    introns = []
    for (chrom, transcript_id, strand), exons in groups.items():
        exons.sort(key=lambda x: x[0])
        gene = exons[0][3] if exons else "."

        for i in range(len(exons) - 1):
            _, _, prev_end = exons[i]
            _, next_start, _ = exons[i + 1]
            intron_start = prev_end
            intron_end = next_start

            if intron_start >= intron_end:
                continue  # Skip overlapping or invalid

            name = f"{gene}:{transcript_id}:intron{i+1}"
            introns.append(
                (chrom, intron_start, intron_end, name, strand, transcript_id)
            )

    # Sort by chrom, start for consistent output
    introns.sort(key=lambda x: (x[0], x[1]))

    # Write output
    out = open(args.output, "w") if args.output else None
    try:
        for chrom, start, end, name, strand, tx_id in introns:
            line = f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n"
            if out:
                out.write(line)
            else:
                print(line, end="")
    finally:
        if out:
            out.close()


if __name__ == "__main__":
    main()
