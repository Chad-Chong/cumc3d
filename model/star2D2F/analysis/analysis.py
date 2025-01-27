############################################
# Plotting utilties for CUMC3D simulations
############################################

#import#
import os
import h5py
from subprocess import call, Popen

################################################################################

print('\t')
print('--------------------------------------------')
print('- Welcome to plotting utilities for CUMC3D -')
print('--------------------------------------------')
print('\t')

################################################################################

# filename#
unsortfile = []

#load#
for root, dirs, files in os.walk('../outfile/'):
  for file in files:
    if file.endswith("nm.hdf5"):
      unsortfile.append(os.path.join(file))

sortedfile = sorted(unsortfile, key=lambda x: int(x.split('-')[1]))
filename = sortedfile[::1]
# filename.append(sortedfile[-1])
print(filename)
# stop now
################################################################################

#arrays for coordinate and dimensions#
n_dim = ['1d', '2d', '3d']
coords = ['cartesian', 'cylindrical', 'spherical']

# read grid parameter files, get model dimension and coordinate #
f = h5py.File('../outfile/'+filename[0], 'r')
dset = f['dimension']
dim = int(dset[:][0])
dset = f['coordinate']
coord = int(dset[:][0])

################################################################################

# assign #
n_dim = n_dim[dim - 1]
coords = coords[coord - 1]

# run process according to user input #
print('\t')
print('Now, run the plotting script ...')
print('\t')
if(n_dim == '1d'):
  script = './plot/plotprofile-1d.py'
else:
  script = './plot/'+str(coords)+'.py'

#run python script, loop over filename #
for i in range (0, len(filename)):
  hdf5_file = '../outfile/'+filename[i]
  Popen(['python3', script, hdf5_file], close_fds=True)

################################################################################

print('-------------------------')
print('- End plotting analysis -')
print('-------------------------')
print('\t')

################################################################################
