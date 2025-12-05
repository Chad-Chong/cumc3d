#!/bin/bash
source /home/cnchong/miniconda3/etc/profile.d/conda.sh
conda activate Rotating_WD
python multiple_compile.py #Get data and compile accordingly CUMC3D
find ./runs -type f -executable -execdir sh -c '{} -> tmp.txt &' \; #Run each CUMC3D and output to tmp.txt