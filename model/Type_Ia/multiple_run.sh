#!/bin/bash
source /home/cnchong/miniconda3/etc/profile.d/conda.sh
conda activate Rotating_WD
python multiple_compile.py #Get data and compile accordingly CUMC3D

count=$(find ./runs -type f -executable | wc -l)
echo "Total number of executables found: $count"

find ./runs -type f -executable -execdir sh -c '{} -> tmp.txt &' \; #Run each CUMC3D and output to tmp.txt