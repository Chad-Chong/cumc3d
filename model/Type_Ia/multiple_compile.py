import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import helmeos
from glob import glob
import subprocess
import os

data_paths = sorted(glob("/home/cnchong/Codes/Helm_RotatingMagnetisedWhiteDwarf/Profile/Output_CentralDensity_*_AxisRatio_*/"))
compile_path = '/home/cnchong/Codes/cumc3d/model/Type_Ia'

for i in data_paths:
    data_path = '/home/cnchong/Codes/cumc3d/model/Type_Ia/runs'+i[61:]
    try: 
        os.mkdir(data_path)
        os.mkdir(data_path+'profile/')
        os.mkdir(data_path+'outfile/')
    except:
        pass
    trim_flag = 0
    for path in os.listdir(i):
        if path == 'trim':
            trim_flag = 1
            break
    if trim_flag == 1:
        subprocess.run(['cp', '-a', i+'trim/.', data_path+'profile/'])
    else:
        subprocess.run(['cp', '-a', i+'.', data_path+'profile/'])
    for file in os.listdir(data_path+'profile'):
        if file == 'hydro_x1_fgrid.dat':
            x1 = np.loadtxt(data_path+'profile/'+file)
            nx = int(x1.shape[0]-7)
        if file == 'hydro_x3_fgrid.dat':
            x3 = np.loadtxt(data_path+'profile/'+file)
            nz = int(x3.shape[0]-7)
    subprocess.run(['sed', '-i', '73s/.*/INTEGER, PARAMETER :: nx = '+str(nx)+'/', compile_path+'/custom-mhd/parameter.h'])
    subprocess.run(['sed', '-i', '75s/.*/INTEGER, PARAMETER :: nz = '+str(nz)+'/', compile_path+'/custom-mhd/parameter.h'])
    subprocess.run(['make', '-j8'])
    subprocess.run(['mv', compile_path+'/CUMC3D', data_path])
    subprocess.run(['make', 'clean'])