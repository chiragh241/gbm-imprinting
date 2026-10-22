#!/usr/bin/env bash
# Create the five conda environments and apply the tool fixes the results depend on.
#
#   gi_analysis     downstream analysis, statistics and figures (Python 3.13, R, bedtools)
#   gi_reads        EM-seq (methylpy, PATCHED below), RNA-seq and ChIP-seq processing, variant
#                   calling, downloads
#   gi_bsseeker2    BS-Seeker2 2.1.7 for the allele-resolved EM-seq, PATCHED (below)
#   gi_cgmaptools   CGmapTools at commit afa6e40 (0.1.3), for SNV and allele-specific calls
#   gi_rice         rice WGBS: Bismark, SNPsplit, and the resequencing SNP calls
#
# The methylpy patch (envs/patches/methylpy.patch) stops methylpy 1.4.7 passing cutadapt options
# that cutadapt 5.2 no longer takes, and gives Picard 3 MarkDuplicates the read group it needs.
# Unpatched, read trimming fails and no library is mapped; see the patch header.
#
# The BS-Seeker2 patch (envs/patches/bs_seeker2-call_methylation.patch) makes the pysam pileup
# keep read pairs without the proper-pair flag and bases without qualities. Unpatched, the tool
# writes empty ATCGmap files and exits 0; see the patch header.
#
# Conda hardlinks environment files to its package cache, so a patch applied inside an
# environment would also change the cached package and every environment linked to it. The
# file is therefore replaced by a copy of itself (a new inode) before it is patched.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"

have() { conda env list | awk '{print $1}' | grep -qx "$1"; }

have gi_analysis || conda env create -q -f "$HERE/gi_analysis.yml"

if ! have gi_reads; then
    conda env create -q -f "$HERE/gi_reads.yml"
    conda activate gi_reads
    site="$(python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
    for f in call_mc_pe.py call_mc_se.py; do
        cp -p "$site/methylpy/$f" "$site/methylpy/$f.copy" && mv -f "$site/methylpy/$f.copy" "$site/methylpy/$f"
    done
    # skipped if the files are already patched (a conda cache patched in place earlier)
    patch -d "$site" -p1 -R -s -f --dry-run < "$HERE/patches/methylpy.patch" > /dev/null 2>&1 \
        || patch -d "$site" -p1 < "$HERE/patches/methylpy.patch"
    conda deactivate
fi

if ! have gi_bsseeker2; then
    conda env create -q -f "$HERE/gi_bsseeker2.yml"
    conda activate gi_bsseeker2
    f="$CONDA_PREFIX/bin/bs_seeker2-call_methylation.py"
    cp -p "$f" "$f.copy" && mv -f "$f.copy" "$f"
    patch -R -s -f --dry-run "$f" "$HERE/patches/bs_seeker2-call_methylation.patch" > /dev/null 2>&1 \
        || patch "$f" "$HERE/patches/bs_seeker2-call_methylation.patch"
    conda deactivate
fi

if ! have gi_cgmaptools; then
    conda env create -q -f "$HERE/gi_cgmaptools.yml"
    conda activate gi_cgmaptools
    src="$CONDA_PREFIX/share/cgmaptools"
    git clone -q https://github.com/guoweilong/cgmaptools.git "$src"
    git -C "$src" checkout -q afa6e4088987da555b4b3482298c93ea0d29c669
    # cgmaptools runs every subcommand from bin/. Its install.sh links the Python, Perl and R
    # tools there and compiles the C tools; only linked tools are used here (snv, asm), so they
    # are linked the same way and nothing is compiled.
    mkdir -p "$src/bin"
    for pr in "$src"/src/*.py "$src"/src/*.pl "$src"/src/*.R; do
        chmod +x "$pr"
        ln -sfn "../src/$(basename "$pr")" "$src/bin/$(basename "$pr" | cut -d. -f1)"
    done
    mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
    echo "export PATH=$src:\$PATH" > "$CONDA_PREFIX/etc/conda/activate.d/cgmaptools.sh"
    conda deactivate
fi

have gi_rice || conda env create -q -f "$HERE/gi_rice.yml"

conda activate gi_reads
grep -q "samtools addreplacerg" \
    "$(python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')/methylpy/call_mc_se.py" \
    || { echo "methylpy patch is not applied" >&2; exit 1; }
conda deactivate
conda activate gi_bsseeker2
grep -q "ignore_orphans=False, min_base_quality=0" \
    "$CONDA_PREFIX/bin/bs_seeker2-call_methylation.py" \
    || { echo "BS-Seeker2 patch is not applied" >&2; exit 1; }
conda deactivate
conda activate gi_cgmaptools
for t in SNVFromATCGmap ASM; do
    [ -x "$(dirname "$(command -v cgmaptools)")/bin/$t" ] \
        || { echo "cgmaptools bin/$t is missing" >&2; exit 1; }
done
conda deactivate
echo "environments ready"
