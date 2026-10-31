#!/bin/sh
# Render the saved chains with the installed CosmoSIS postprocessor.
set -eu

project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
postprocess="${COSMOSIS_POSTPROCESS:-cosmosis-postprocess}"
output_dir="$project_dir/output/curved/native_postprocess"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/cosmosis_native_figures.XXXXXX")
trap 'rm -rf "$scratch"' EXIT HUP INT TERM
export MPLCONFIGDIR="$scratch/mpl"
export PYTHONDONTWRITEBYTECODE=1

render() {
    group="$1"
    destination="$output_dir/$group"
    legend="$2"
    shift 2
    mkdir -p "$destination"
    for format in png; do
        if ! COSMOSIS_FIGURE_GROUP="$group" "$postprocess" --burn 0.3 --no-2d \
            --file-type "$format" --tweaks "$project_dir/scripts/native_plot_tweaks.py" \
            --outdir "$scratch/postprocess" --legend "$legend" "$@" \
            > "$scratch/postprocess.log" 2>&1; then
            cat "$scratch/postprocess.log" >&2
            exit 1
        fi
        cp "$scratch/postprocess/cosmological_parameters--h0.$format" "$destination/h0.$format"
        cp "$scratch/postprocess/cosmological_parameters--w.$format" "$destination/w.$format"
        rm -rf "$scratch/postprocess"
    done
    printf 'Created %s\n' "$destination"
}

cd "$project_dir"
chains=.run-data/curved
if [ "${1:-}" != "--corner-only" ]; then
render omega_k_minus_0p01 'Pantheon+SH0ES (Omega_k=-0.01)|DES-Y3 BAO (Omega_k=-0.01)' \
    "$chains/pantheonplussh0es_ok_-0.01.txt" "$chains/des_y3_bao_ok_-0.01.txt"
render omega_k_zero 'Pantheon+SH0ES (Omega_k=0.00)|DES-Y3 BAO (Omega_k=0.00)' \
    "$chains/pantheonplussh0es_ok_+0.00.txt" "$chains/des_y3_bao_ok_+0.00.txt"
render omega_k_plus_0p01 'Pantheon+SH0ES (Omega_k=+0.01)|DES-Y3 BAO (Omega_k=+0.01)' \
    "$chains/pantheonplussh0es_ok_+0.01.txt" "$chains/des_y3_bao_ok_+0.01.txt"
render pantheon_curvature 'Omega_k=-0.01|Omega_k=0.00|Omega_k=+0.01' \
    "$chains/pantheonplussh0es_ok_-0.01.txt" "$chains/pantheonplussh0es_ok_+0.00.txt" \
    "$chains/pantheonplussh0es_ok_+0.01.txt"
render bao_curvature 'Omega_k=-0.01|Omega_k=0.00|Omega_k=+0.01' \
    "$chains/des_y3_bao_ok_-0.01.txt" "$chains/des_y3_bao_ok_+0.00.txt" \
    "$chains/des_y3_bao_ok_+0.01.txt"
fi

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
