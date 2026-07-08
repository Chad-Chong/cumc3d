#!/bin/bash
source /home/cnchong/miniconda3/etc/profile.d/conda.sh
conda activate Rotating_WD
python multiple_compile.py #Get data and compile accordingly CUMC3D

SEARCH_PATH="/mnt3/cnchong/thesis_runs"

count=$(find "$SEARCH_PATH" -type f -name "CUMC3D" -executable | wc -l)
echo "Total number of executables found: $count"

find "$SEARCH_PATH" -type f -name "CUMC3D" -execdir sh -c './CUMC3D > tmp.txt 2>&1 &' \;