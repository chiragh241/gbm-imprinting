# Paths and settings shared by every script. Source it; do not run it.
#
#   GI_RUN      one directory for everything a run downloads or writes (default ~/gbm_imprinting_run)
#   GI_DATA     downloaded inputs: genomes, annotations, published tables, reads
#   GI_WORK     intermediate files
#   GI_RESULTS  final tables and figures
#   GI_THREADS  threads for aligners and sorters
#
# Any of them can be set in the environment before sourcing to override the default.
GI_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GI_RUN="${GI_RUN:-$HOME/gbm_imprinting_run}"
GI_DATA="${GI_DATA:-$GI_RUN/data}"
GI_WORK="${GI_WORK:-$GI_RUN/work}"
GI_RESULTS="${GI_RESULTS:-$GI_RUN/results}"
GI_THREADS="${GI_THREADS:-8}"
GI_SHIPPED="$GI_REPO/data"
export GI_REPO GI_RUN GI_DATA GI_WORK GI_RESULTS GI_THREADS GI_SHIPPED

# Activate one of the conda environments defined in envs/.
gi_env() {
    if [ -z "${CONDA_EXE:-}" ]; then
        echo "conda is not on PATH; install it and create the environments in envs/" >&2
        return 1
    fi
    # Activation scripts read unset variables, so nounset is lifted around them.
    local had_u=0
    case $- in *u*) had_u=1; set +u ;; esac
    # shellcheck disable=SC1091
    source "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
    conda activate "$1"
    local rc=$?
    [ "$had_u" = 1 ] && set -u
    return $rc
}

# Download a URL to a path unless the path already exists, then check it is not empty.
gi_fetch() {
    local url="$1" out="$2"
    if [ -s "$out" ]; then return 0; fi
    mkdir -p "$(dirname "$out")"
    wget -q -O "$out.part" "$url" && mv -f "$out.part" "$out"
    [ -s "$out" ] || { echo "download failed: $url" >&2; return 1; }
}
