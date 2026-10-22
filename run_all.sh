#!/usr/bin/env bash
# Run the pipeline in order: every numbered script, module by module.
#
#   ./run_all.sh                                   everything
#   ./run_all.sh scripts/02_gene_sets              one module
#   ./run_all.sh scripts/02_gene_sets/03_build_sets.py   one script
#
# Python and R scripts run in the gi_analysis environment; shell scripts activate the
# environment they need themselves. Each script's output goes to $GI_RUN/logs.
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$GI_RUN/logs"

run_one() {
    local s="$1"
    local log="$GI_RUN/logs/$(basename "$(dirname "$s")")__$(basename "$s").log"
    echo "=== $(date '+%F %T') $s"
    case "$s" in
        *.py) ( gi_env gi_analysis && python "$s" ) ;;
        *.R)  ( gi_env gi_analysis && Rscript "$s" ) ;;
        *.sh) bash "$s" ;;
        *) echo "unknown script type: $s" >&2; return 1 ;;
    esac > "$log" 2>&1 || { echo "failed: $s (log: $log)" >&2; return 1; }
}

targets=("${@:-$GI_REPO/scripts}")
for t in "${targets[@]}"; do
    if [ -d "$t" ]; then
        while IFS= read -r s; do run_one "$s"; done < <(
            find "$t" -type f -regex '.*/[0-9][0-9]_[^/]*\.\(py\|R\|sh\)' | sort)
    else
        run_one "$t"
    fi
done
echo "=== $(date '+%F %T') done"
