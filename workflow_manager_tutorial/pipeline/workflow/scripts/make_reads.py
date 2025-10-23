#!/usr/bin/env python3
import sys
from pathlib import Path
import random

out = Path(snakemake.output.fastq)
out.parent.mkdir(parents=True, exist_ok=True)

# read reference and sample id
ref = snakemake.input[0]
sample = snakemake.wildcards.sample
n = int(snakemake.params.reads_per_sample)

seqs = []
with open(ref) as fh:
        current = ""
        for line in fh:
                if line.startswith(">"):
                        if current:
                                seqs.append(current)
                                current = ""
                else:
                        current += line.strip()
        if current:
                seqs.append(current)

# generate n single-end reads by sampling random substrings (length 50)
with open(out, "w") as fh:
        for i in range(n):
                s = random.choice(seqs)
                if len(s) < 60:
                        pos = 0
                        read = s
                else:
                        pos = random.randint(0, len(s) - 50)
                        read = s[pos:pos+50]
                fh.write(f"@{sample}_{i}\\n{read}\\n+\\n{'I'*len(read)}\\n")
