#!/bin/bash

# Define paths
previous_run_directory="/mnt2/shcheung/MagnetarFormationHydro"
parent_directory="/mnt2/shcheung/MagnetarFormationHydro"
new_run_directory_prefix="LB2_CentralDensity_10.699_AxisRatio"
template_path="/home/shcheung/AIC_problems/AIC_template"

# Create the parent directory if it doesn't exist
mkdir -p "$parent_directory"

# Loop through each previous run directory
for old_run_directory in "$previous_run_directory"/LB_CentralDensity_10.699_AxisRatio_*_alpha_2d-9_kappa_1d-11; do
    # Extract the AxisRatio from the directory name
    axis_ratio=$(basename "$old_run_directory" | sed -n 's/.*AxisRatio_\([0-9.]*\)_alpha_2d-9_kappa_1d-11/\1/p')
    
    # Define the new hydro run directory
    new_hydro_directory="$parent_directory/${new_run_directory_prefix}_${axis_ratio}_alpha_2d-9_kappa_1d-11"
    
    # Create the new hydro directory
    mkdir -p "$new_hydro_directory"
    
    # Copy the template files to the new hydro directory
    cp -r "$template_path"/* "$new_hydro_directory/"

    # Create the profile directory in the new hydro directory
    mkdir -p "$new_hydro_directory/profile/Profile"
    
    # Copy the .hdf5 file from the previous run directory to the new hydro directory
    cp "$old_run_directory/outfile/rkiter-65-nm.hdf5" "$new_hydro_directory/profile/"
    cp "$old_run_directory/profile/hydro_x1_fgrid.dat" "$new_hydro_directory/profile/"
    cp "$old_run_directory/profile/hydro_x2_fgrid.dat" "$new_hydro_directory/profile/"

    # Navigate to the new hydro directory and start the simulation
    pushd "$new_hydro_directory" > /dev/null
    make -j8
    nohup ./CUMC3D > tmp.txt &
    popd > /dev/null
done