#!/usr/bin/env python3
import sys
import pysam
from pathlib import Path

bam = snakemake.input.bam
out = Path(snakemake.output[0])
out.parent.mkdir(parents=True, exist_ok=True)

samfile = pysam.AlignmentFile(bam, "rb")
total = 0
mapped = 0
for r in samfile.fetch(until_eof=True):
        total += 1
        if not r.is_unmapped:
                mapped += 1
samfile.close()

with open(out, "w") as fh:
        fh.write(f"total_reads\\t{total}\\n")
        fh.write(f"mapped_reads\\t{mapped}\\n")
        fh.write(f"mapping_rate\\t{mapped/total if total else 0:.3f}\\n")
