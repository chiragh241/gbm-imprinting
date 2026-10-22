"""Write a copy of a genome with every SNP position set to N.

usage: python nmask.py genome.fa snps.txt out.fa

snps.txt is the SNPsplit SNP file from 03_call_snps.sh (ID, Chr, Position (1-based), SNP value,
Ref/SNP). The reference base at each position is checked against the file's Ref before masking,
so a coordinate slip would stop the script rather than mask the wrong bases.
"""
import sys
from collections import defaultdict

genome, snps, out = sys.argv[1:4]

pos = defaultdict(list)
ref = {}
with open(snps) as fh:
    next(fh)
    for line in fh:
        _, chrom, p, _, alleles = line.rstrip("\n").split("\t")
        p = int(p) - 1
        pos[chrom].append(p)
        ref[(chrom, p)] = alleles.split("/")[0].upper()


def emit(name, seq, fo):
    n = 0
    for p in pos.get(name, []):
        if chr(seq[p]).upper() != ref[(name, p)]:
            sys.exit(f"{name}:{p + 1} is {chr(seq[p])} in the genome, {ref[(name, p)]} in the SNP file")
        seq[p] = ord("N")
        n += 1
    fo.write(f">{name}\n")
    for i in range(0, len(seq), 60):
        fo.write(seq[i:i + 60].decode() + "\n")
    print(f"{name}\t{len(seq)}\t{n} masked")
    return n


total = 0
seen = set()
with open(genome) as fi, open(out, "w") as fo:
    name, parts = None, []
    for line in fi:
        if line.startswith(">"):
            if name is not None:
                total += emit(name, bytearray(b"".join(parts)), fo)
                seen.add(name)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip().encode())
    total += emit(name, bytearray(b"".join(parts)), fo)
    seen.add(name)
missing = set(pos) - seen
if missing:
    sys.exit(f"SNPs on sequences not in the genome: {sorted(missing)}")
print(f"total\t{total} masked of {sum(len(v) for v in pos.values())} SNPs")
