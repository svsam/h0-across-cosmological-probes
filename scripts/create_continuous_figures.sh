#!/bin/sh
# Render only the GetDist corner plots from the saved continuous chains.
set -eu
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
"${COSMOSIS_PYTHON:-python}" "$project_dir/scripts/getdist_corner.py"
