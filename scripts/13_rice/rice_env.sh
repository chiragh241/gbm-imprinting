# Paths and settings shared by the rice shell scripts. Source it; do not run it.
source "$(dirname "${BASH_SOURCE[0]}")/../../config.sh"

RICE=$GI_DATA/rice
REF=$RICE/ref                   # IRGSP-1.0 with ChrC and ChrM, MSU7 gene models (01)
SUPP=$RICE/supplementary        # Tonosaki et al. 2024 Supplementary Data 1 and 3 (01)
SRA=$RICE/sra                   # raw runs, .sra (02)
RW=$GI_WORK/13_rice             # everything derived
RUNS=$GI_SHIPPED/rice/runs.tsv  # run, type (wgbs or reseq), alt parent, mother, father, replicate

# Large temporary files (sorts, Bismark) stay off /tmp, which is often a small tmpfs.
mkdir -p "$RW/tmp"
export TMPDIR=$RW/tmp

# For testing only: with GI_TEST_SPOTS=N every run is read up to its first N spots (read pairs).
SPOTS=()
[ -n "${GI_TEST_SPOTS:-}" ] && SPOTS=(-X "$GI_TEST_SPOTS")

# Bismark alignment options. Tonosaki et al. allowed "1 mismatch within 32-nt regions in seeds"
# (bowtie2 -N 1 -L 32).
BISMARK_OPTS="-N 1 -L 32"

# Trim Galore's defaults (quality 20, Illumina adapter on both reads, 1 bp overlap, 10% error,
# pairs dropped when either read is under 20 bp), run as one paired cutadapt call so the reads can
# stream from the .sra without an uncompressed FASTQ on disk.
TRIM_ARGS="-q 20 -O 1 -e 0.1 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 20"

# Phred offset of a run's qualities: 64 when no quality character below '@' appears in its first
# 20,000 spots, else 33. The resequencing runs are Phred+64 and would otherwise go untrimmed; the
# WGBS runs are Phred+33. The lowest character is computed over the whole sample: an early-exiting
# grep here would raise SIGPIPE upstream, which pipefail turns into a false "64".
qual_base() {
    local low
    low=$(fastq-dump --stdout -X 20000 "$SRA/$1.sra" 2>/dev/null | awk '
        BEGIN { for (i = 33; i < 127; i++) ord[sprintf("%c", i)] = i; m = 999 }
        NR % 4 == 0 { for (j = 1; j <= length($0); j++) { c = ord[substr($0, j, 1)]; if (c < m) m = c } }
        END { print m }')
    if [ -z "$low" ] || [ "$low" = 999 ]; then echo "qual_base: no reads from $1" >&2; return 1; fi
    if [ "$low" -ge 64 ]; then echo 64; else echo 33; fi
}
