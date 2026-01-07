#!/bin/bash
read -p "Enter directory: " directory

profile=${directory%%/}/profile
outfile=${directory%%/}/outfile

mv CUMC3D $directory
mkdir -p $profile
cp ./profile/get_test.sh $profile
mkdir -p $outfile