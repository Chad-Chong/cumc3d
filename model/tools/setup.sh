#!/bin/bash

# Define paths
initial_data_path="/home/shcheung/RotatingMagnetisedWD/git_RMWD/Profile"
template_path="/home/shcheung/AIC_problems/AIC_template"
# source_path="/home/shcheung/AIC_problems/model/AIC_src"
parent_directory="/mnt2/shcheung/MagnetarFormationHydro"
python_script="/home/shcheung/RotatingMagnetisedWD/git_RMWD/Convert2Hydro.py"  # Update this path to the actual location of Convert2Hydro.py

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate star

# Create the parent directory if it doesn't exist
mkdir -p "$parent_directory"

# cp -r "$source_path"/ "$parent_directory/"
# Loop through each set of initial data files
for grid_file in "$initial_data_path"/Output_CentralDensity_10.699_AxisRatio_*_grid.dat; do
    # Extract the AxisRatio from the filename
    axis_ratio=$(basename "$grid_file" | sed -n 's/.*AxisRatio_\([0-9.]*\)_grid.dat/\1/p')
    
    # Define the corresponding hydro file
    hydro_file="$initial_data_path/Output_CentralDensity_10.699_AxisRatio_${axis_ratio}_hydro.dat"
    
    # Define the directory name for the hydro simulation
    hydro_directory="$parent_directory/LB_CentralDensity_10.699_AxisRatio_${axis_ratio}_alpha_2d-9_kappa_1d-11"
    
    # Create the hydro directory
    mkdir -p "$hydro_directory"
    
    # Copy the template files to the hydro directory
    cp -r "$template_path"/* "$hydro_directory/"

    mkdir -p "$hydro_directory/profile/Profile"
    # Copy the initial data files to the hydro directory
    cp "$grid_file" "$hydro_directory/profile/Profile/"
    cp "$hydro_file" "$hydro_directory/profile/Profile/"
    
    mkdir -p "$hydro_directory/profile/Plot"
    pushd "$hydro_directory/profile" > /dev/null
    # Run the Python script to generate the hydro .dat files
    python3 "$python_script" "$hydro_directory/profile/Output_CentralDensity_10.699_AxisRatio_${axis_ratio}_grid.dat" "$hydro_directory/profile/Output_CentralDensity_10.699_AxisRatio_${axis_ratio}_hydro.dat"
    
    # Move the generated .dat files to the profile directory
    # mv hydro_*.dat "$hydro_directory/profile/"

    pushd "$hydro_directory" > /dev/null
    make -j8
    # nohup ./CUMC3D > tmp.txt &
    popd > /dev/null
    
    # Return to the previous directory
    popd > /dev/null
done

echo "All hydro simulations have been set up."
