import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import helmeos
from glob import glob
import subprocess
import os
import math as m

#####################################################################
#####################################################################

gconst = 6.67430e-8
clight = 2.99792458e10
solar = 1.98847e33
lencgs2code = (clight**2)/(solar*gconst)
masscgs2code = (1.0e0/solar)
rhocgs2code = (masscgs2code/lencgs2code**3)
tcgs2code = (clight**3)/(solar*gconst)
energycgs2code = (1.0E0/clight**2)
avo = 6.0221367e23
amu = 1.6605402e-24
ev2erg = 1.60217648740e-12
mev2erg = ev2erg*1.0e6
mev2gr = mev2erg/clight**2

def filtering(array,atmo_indicies):
    dummy = np.copy(array)
    dummy[atmo_indicies] = 0
    return dummy

def helmeos_e(dens, temp, abar, zbar):
    out = helmeos.eos_DT(dens, temp, abar, zbar)
    return out['etot']

# set the number of nucleons in the element
aion = np.array([4,12,16,20,14,28,56])

# set the number of protons in the element
zion = np.array([2,6,8,10,12,14,28])

# set the binding energy of the element
bion = np.array([28.29603,92.16294,127.62093,160.64788,198.25790,236.53790,484.00300])

# ! set the number of neutrons and mass
nion = aion - zion


# mass of each isotope
mion = nion*1.67492721184e-24 + zion*1.67262163783e-24 - bion*mev2gr          

# a common approximation   
wion = aion       

def azbar(xmass):
    # molar abundances
    ymass = xmass.T/wion

    # mean molar mass
    wbar  = 1.0e0/np.sum(ymass, axis = 2)

    # mean number of nucleons
    sum1  = np.sum(aion*ymass, axis=2)
    abar  = wbar * sum1

    # mean charge
    ye  = np.sum(zion*ymass, axis=2)
    zbar  = wbar * ye

    return abar.T, zbar.T

def azbar_s(xmass):
    # molar abundances
    ymass = xmass/wion

    # mean molar mass
    wbar  = 1.0e0/np.sum(ymass)

    # mean number of nucleons
    sum1  = np.sum(aion*ymass)
    abar  = wbar * sum1

    # mean charge
    ye  = np.sum(zion*ymass)
    zbar  = wbar * ye

    return abar, zbar

irho = 0
ivx = 1
ivy = 2
ivz = 3
itau = 4
ihe4 = 5
ic12 = 6
io16 = 7
ine20 = 8
img24 = 9
isi28 = 10
ini56 = 11
iye2 = 12
iturbq = 13
iscaG1 = 14
iscaG2 = 15

#####################################################################
#####################################################################

working_directory = './outfile'
profile_directory = './profile'
print('Enter the index of the file for mapping to SNEC:')
n = input()

# Read file
filename = working_directory+"/rkiter-"+str(0)+"-nm.hdf5"
f0 = h5py.File(filename, "r")
primitive = f0['primitive'][:].T
atmo = f0['prim_a'][irho]
rho = primitive[0,1:-1,1,1:-1]
vol = (f0['vol'][:].T)[:,0,:]
atmo_indicies = rho>atmo
initial_mass = np.sum((rho*vol)[atmo_indicies])

# Face value
x1f = np.loadtxt(profile_directory+'/hydro_x1_fgrid.dat').T
x3f = np.loadtxt(profile_directory+'/hydro_x3_fgrid.dat').T
# Center value
x1f = x1f [3:-3]*lencgs2code # Remove BC
x1c = (x1f[1:]+x1f[:-1])/2 #mid pt
x3f = x3f[3:-3]*lencgs2code # Remove BC
x3c = (x3f[1:]+x3f[:-1])/2 #mid pt

xx,yy= np.meshgrid(x1c,x3c,indexing='ij')
rr_grid = np.sqrt(xx**2+yy**2)

xxf, yyf = np.meshgrid(x1f,x3f,indexing='ij')

nr = m.floor((len(x1f)-1)/50 -1)*50
dr = (x1f[-1]-x1f[0])/nr
dx = dr/2 #r coordinates in the center
vol = (f0['vol'][:].T)[:,0,:]
r_arr = np.linspace(dx, x1f[-1]-dx, num=nr)

filename = working_directory+"/rkiter-"+str(n)+"-nm.hdf5"
f = h5py.File(filename, "r")

primitive = f['primitive'][:].T
eps = (f['epsilon'][:].T)[1:-1,1,1:-1]
prim_a = f['prim_a'][:]
atmo = f['prim_a'][irho]
temp = (f['temp'][:].T)[1:-1,1,1:-1]
abar,zbar = azbar(primitive[ihe4:ini56+1, 1:-1,1,1:-1])

rho = primitive[irho,1:-1,1,1:-1]
atmo_indicies = rho<=atmo
v1 = primitive[ivx,1:-1,1,1:-1]
v2 = primitive[ivy,1:-1,1,1:-1]
v3 = primitive[ivz,1:-1,1,1:-1]
Ye = primitive[iye2,1:-1,1,1:-1]

Den_he = primitive[ihe4,1:-1,1,1:-1]*rho
Den_c = primitive[ic12,1:-1,1,1:-1]*rho
Den_o = primitive[io16,1:-1,1,1:-1]*rho
Den_ne = primitive[ine20,1:-1,1,1:-1]*rho
Den_mg = primitive[img24,1:-1,1,1:-1]*rho
Den_si = primitive[isi28,1:-1,1,1:-1]*rho
Den_ni = primitive[ini56,1:-1,1,1:-1]*rho

#Number density of nucleons / baryon and electron
Num_he = primitive[ihe4,1:-1,1,1:-1]*rho/rhocgs2code/mion[0]*avo*zion[0]
Num_c = primitive[ic12,1:-1,1,1:-1]*rho/rhocgs2code/mion[1]*avo*zion[1]
Num_o = primitive[io16,1:-1,1,1:-1]*rho/rhocgs2code/mion[2]*avo*zion[2]
Num_ne = primitive[ine20,1:-1,1,1:-1]*rho/rhocgs2code/mion[3]*avo*zion[3]
Num_mg = primitive[img24,1:-1,1,1:-1]*rho/rhocgs2code/mion[4]*avo*zion[4]
Num_si = primitive[isi28,1:-1,1,1:-1]*rho/rhocgs2code/mion[5]*avo*zion[5]
Num_ni = primitive[ini56,1:-1,1,1:-1]*rho/rhocgs2code/mion[6]*avo*zion[6]
Num_n = Num_he + Num_c + Num_o + Num_ne + Num_mg + Num_si + Num_ni
Num_e = Num_n * (np.minimum(Ye*abar, zbar)/abar)

#####################################################################
#####################################################################

masshell_r = []
mass_r = []
ke_r = []
dens_r = []
vel_r = []
inte_r = []
abar_r = []
zbar_r = []
Ye_r = []
Tguess_r = []
xiso_r = []
r_vol_r = []
ir = -1

for r in r_arr:
    ir += 1
    filtered_rho = filtering(rho, atmo_indicies)
    r_indicies = np.logical_and(r-dx<=rr_grid, rr_grid<r+dx)
    inner_r_indicies = (rr_grid <= (r-dx))
    r_vol = 4*np.pi*((r+dr)**3-(r-dr)**3)/3

    mass_coord = np.sum((filtered_rho*vol)[inner_r_indicies])
    mass = np.sum((filtered_rho*vol)[r_indicies])
    if mass == 0:
        flag = True #Reached a radius where there is only atmosphere, no need to look further
        break
    ke = 0.5 * np.sum( (filtered_rho*vol*(v1**2+v3**2) )[r_indicies] )
    inte = np.sum((filtered_rho*vol*eps )[r_indicies] )/mass #energy per mass in this shell
    dens = mass/r_vol

    # Xiso
    filtered = filtering(Den_he, atmo_indicies)
    xiso_he = np.sum((filtered*vol)[r_indicies])/r_vol/dens

    filtered = filtering(Den_c, atmo_indicies)
    xiso_c = np.sum((filtered*vol)[r_indicies])/r_vol/dens
    
    filtered = filtering(Den_o, atmo_indicies)
    xiso_o = np.sum((filtered*vol)[r_indicies])/r_vol/dens

    filtered = filtering(Den_ne, atmo_indicies)
    xiso_ne = np.sum((filtered*vol)[r_indicies])/r_vol/dens

    filtered = filtering(Den_mg, atmo_indicies)
    xiso_mg = np.sum((filtered*vol)[r_indicies])/r_vol/dens

    filtered = filtering(Den_si, atmo_indicies)
    xiso_si = np.sum((filtered*vol)[r_indicies])/r_vol/dens

    filtered = filtering(Den_ni, atmo_indicies)
    xiso_ni = np.sum((filtered*vol)[r_indicies])/r_vol/dens

    xiso = np.array([xiso_he, xiso_c, xiso_o, xiso_ne, xiso_mg, xiso_si, xiso_ni])

    # Ye
    filtered = filtering(Num_e, atmo_indicies)
    Ne_r = np.sum((filtered*vol)[r_indicies])

    filtered = filtering(Num_n, atmo_indicies)
    Ne_n = np.sum((filtered*vol)[r_indicies])

    Ye = Ne_r/Ne_n

    # Abar and Zbar

    abar, zbar = azbar_s(xiso)

    zbar = np.minimum(zbar, Ye*abar)

    # Guess Temperature

    min_index =  np.argmin(np.abs((filtered_rho-dens)[r_indicies]))

    masshell_r.append(mass/masscgs2code)
    mass_r.append(mass_coord/masscgs2code)
    dens_r.append(dens/rhocgs2code)
    vel_r.append(np.sqrt(2*ke/mass)*clight)
    ke_r.append(ke/masscgs2code*clight**2)
    inte_r.append(inte/energycgs2code)
    abar_r.append(abar)
    zbar_r.append(zbar)
    Ye_r.append(Ye)
    xiso_r.append(xiso)
    Tguess_r.append(temp[r_indicies][min_index]*1e9)
    r_vol_r.append(r_vol)

masshell_r = np.array(masshell_r)
mass_r = np.array(mass_r)
ke_r = np.array(ke_r)
dens_r = np.array(dens_r)
vel_r = np.array(vel_r)
inte_r = np.array(inte_r)
abar_r = np.array(abar_r)
zbar_r = np.array(zbar_r)
Ye_r = np.array(Ye_r)
Tguess_r = np.array(Tguess_r)
xiso_r = np.array(xiso_r).T
r_vol_r = np.array(r_vol_r)

all = np.array((inte_r,dens_r,abar_r,zbar_r,Ye_r,Tguess_r))
r_arr = r_arr[:ir]

#####################################################################
#####################################################################

output_n = np.zeros(1, dtype='int')
output_n[0] = len(r_arr)
np.savetxt('./Inter_SNEC.dat', output_n, fmt='%i')
f = open('./Inter_SNEC.dat', 'a')
np.savetxt(f, all)
f.close()

#####################################################################
#####################################################################

subprocess.run('./Invert_Helm.out')

#####################################################################
#####################################################################

T_r = np.loadtxt('./Python_SNEC.dat')

#####################################################################
#####################################################################

print("Enter the file name for SNEC")
file_name = input()
indicies = np.arange(output_n[0])+1
omega = np.zeros(output_n[0])
SNEC_arr = np.array((indicies, mass_r, r_arr/lencgs2code, T_r, dens_r, vel_r, Ye_r, omega))
fmt_arr = ['%i', '%1.18e', '%1.18e', '%1.18e', '%1.18e', '%1.18e', '%1.18e', '%1.18e']
np.savetxt(working_directory+'/'+file_name+'.short', output_n, fmt='%i')
f = open(working_directory+'/'+file_name+'.short', 'a')
np.savetxt(f, SNEC_arr.T,fmt=fmt_arr, delimiter='     ')
f.close()
os.remove('./Inter_SNEC.dat')
os.remove('./Python_SNEC.dat')

f = open(working_directory+'/'+file_name+'.iso.dat', 'w')
comp_profile_n = np.array((output_n[0], 7))

for i in comp_profile_n:
    f.write(str(i)+'     ')
f.write("\n")


for i in aion:
    f.write(str(i)+'     ')
f.write("\n")

for i in zion:
    f.write(str(i)+'     ')
f.write("\n")

SNEC_profile_arr = np.array((mass_r, r_arr/lencgs2code, xiso_r[0],xiso_r[1],xiso_r[2],xiso_r[3],xiso_r[4],xiso_r[5],xiso_r[6]))
np.savetxt(f, SNEC_profile_arr.T, delimiter='     ')
f.close()

#####################################################################
#####################################################################

print("Number of grids generated are", str(nr))
if flag:
    print("Encountered the atmosphere at index", ir)
print('Initial mass is',initial_mass, 'solar mass')
print('Mass of Ni is',np.sum((xiso_r[6,:]*dens_r*rhocgs2code*r_vol_r)), 'solar mass')
print('Mass coordinate of boundary Ni is', mass_r[-1]*masscgs2code, 'solar mass')