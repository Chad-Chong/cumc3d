######################################################
#
# make.inc file incluing compilation flags for CUMC3D
#
######################################################

# Want to fix singularies pole? #
FIXPOLE = -DFIXPOLE 

# Want to debug? #
#DEBUG = -DDEBUG

# Force legacy behaviour, turn this off when bugged #
#LEGACY = -fpp -vms -nofor-main

# Want to get divergence of B? #
#DIVB = -DDIVB

# Want to look for compilation message? #
#COMPILE = -Minfo -Mipa=fast

# Choose between CPU or GPU threading #
#MULTICORE = -DGPU -acc
MULTICORE = -qopenmp

# Define F90 compiler #
# dummy := $(shell module load hdf5/1.14.0 2>/dev/null)
F90 := h5pfc

# need extra math packages? #
#MATHPACK = -llapack -lblas

# Compiler flags
# PGI = -fastsse -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -O4
# PGI = -g -O2 -Mbounds -Mchkptr -Mchkstk -traceback # for debug
# ifortflags = -g -O3 -traceback -warn all -check all
#ifortomp = -qopenmp
# ifortdebug = -check all -fpe0

# ifort = $(ifortflags) $(ifortomp) $(ifortdebug)

# Define F90 flags #
# -Mipa=fast is suppressed #
F90FLAGS= -cpp -module . $(LEGACY) $(PGI) $(ifort) $(MULTICORE) $(COMPILE) $(DEBUG) $(DIVB) $(FIXPOLE) $(MATHPACK)
LDFLAGS= -cpp -module . $(LEGACY) $(PGI) $(ifort) $(MULTICORE) $(COMPILE) $(DEBUG) $(DIVB) $(FIXPOLE) $(MATHPACK)
