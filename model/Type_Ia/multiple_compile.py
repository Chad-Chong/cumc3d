import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import helmeos
from glob import glob
import subprocess
import os

data_paths = sorted(glob("/mnt3/cnchong/thesis_runs/Output_CentralDensity_*_AxisRatio_*/"))
compile_path = '/home/cnchong/Codes/cumc3d/model/Type_Ia'

for i in data_paths:
    try: 
        subprocess.run(['mkdir', i+'outfile'])
        print("Made", i+"outfile")
    except:
        pass
    for file in os.listdir(i+'profile'):
        if file == 'hydro_x1_fgrid.dat':
            x1 = np.loadtxt(i+'profile/'+file)
            nx = int(x1.shape[0]-7)
        if file == 'hydro_x3_fgrid.dat':
            x3 = np.loadtxt(i+'profile/'+file)
            nz = int(x3.shape[0]-7)
    subprocess.run(['make', 'clean'])
    subprocess.run(['sed', '-i', '73s/.*/INTEGER, PARAMETER :: nx = '+str(nx)+'/', compile_path+'/custom-mhd/parameter.h'])
    subprocess.run(['sed', '-i', '75s/.*/INTEGER, PARAMETER :: nz = '+str(nz)+'/', compile_path+'/custom-mhd/parameter.h'])
    subprocess.run(['make', 'clean'])
    subprocess.run(['make', '-j8'])
    subprocess.run(['mv', compile_path+'/CUMC3D', i])
    subprocess.run(['make', 'clean'])