#!/bin/sh
set -eu
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
postprocess="${COSMOSIS_POSTPROCESS:-cosmosis-postprocess}"
destination="$project_dir/output/continuous/figures"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/cosmosis_continuous_figures.XXXXXX")
trap 'rm -rf "$scratch"' EXIT HUP INT TERM
export MPLCONFIGDIR="$scratch/mpl"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$destination"
for format in png; do
    "$postprocess" --burn 0.3 --thin 10 --no-plots --file-type "$format" \
        --extra "$project_dir/scripts/continuous_postprocess.py" \
        --outdir "$scratch/postprocess" \
        "$project_dir/.run-data/continuous/des_y3_bao.txt" \
        "$project_dir/.run-data/continuous/pantheonplussh0es.txt" \
        "$project_dir/.run-data/continuous/desi_dr2_bao.txt" \
        "$project_dir/.run-data/continuous/planck_compressed.txt"
    cp "$scratch/postprocess/continuous_corner.$format" "$destination/continuous_corner.$format"
    cp "$scratch/postprocess/continuous_h0.$format" "$destination/continuous_h0.$format"
    cp "$scratch/postprocess/continuous_h0_histograms.$format" "$destination/continuous_h0_histograms.$format"
    cp "$scratch/postprocess/continuous_h0_histograms_50_80.$format" "$destination/continuous_h0_histograms_50_80.$format"
    cp "$scratch/postprocess/continuous_parameter_histograms.$format" "$destination/continuous_parameter_histograms.$format"
    rm -rf "$scratch/postprocess"
done
