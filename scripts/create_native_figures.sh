#!/bin/sh
# Render the fixed-curvature corner plot from the saved chains.
set -eu

project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
postprocess="${COSMOSIS_POSTPROCESS:-cosmosis-postprocess}"
output_dir="$project_dir/output/curved/native_postprocess"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/cosmosis_native_figures.XXXXXX")
trap 'rm -rf "$scratch"' EXIT HUP INT TERM
export MPLCONFIGDIR="$scratch/mpl"
export PYTHONDONTWRITEBYTECODE=1

cd "$project_dir"
chains=.run-data/curved
# CosmoSIS runs this supported extension on each of the six chain inputs.
# It draws one corner figure with a panel for each dataset and colors for the
# three fixed curvature choices. A light thinning keeps KDE plotting quick.
destination="$output_dir/corner"
mkdir -p "$destination"
for format in png; do
    if ! "$postprocess" --burn 0.3 --thin 10 --no-plots --file-type "$format" \
        --extra "$project_dir/scripts/corner_postprocess.py" \
        --outdir "$scratch/postprocess" \
        "$chains/pantheonplussh0es_ok_-0.01.txt" \
        "$chains/pantheonplussh0es_ok_+0.00.txt" \
        "$chains/pantheonplussh0es_ok_+0.01.txt" \
        "$chains/des_y3_bao_ok_-0.01.txt" \
        "$chains/des_y3_bao_ok_+0.00.txt" \
        "$chains/des_y3_bao_ok_+0.01.txt" \
        > "$scratch/postprocess.log" 2>&1; then
        cat "$scratch/postprocess.log" >&2
        exit 1
    fi
    cp "$scratch/postprocess/h0_w_corner.$format" "$destination/h0_w_corner.$format"
    rm -rf "$scratch/postprocess"
done
printf 'Created %s\n' "$destination"
