######################################################
#
# make.inc file incluing compilation flags for CUMC3D
#
######################################################

# Want to fix singularies pole? #
FIXPOLE = -DFIXPOLE 

# Want to debug? #
#DEBUG = -DDEBUG

# Want to get divergence of B? #
#DIVB = -DDIVB

# Want to look for compilation message? #
#COMPILE = -Minfo -Mipa=fast

# Choose between CPU or GPU threading #
#MULTICORE = -DGPU -acc
#MULTICORE = -mp

# Define F90 compiler #
F90=h5pfc

# need extra math packages? #
#MATHPACK = -llapack -lblas

# PGI Compiler flags
#PGI = -fastsse -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -O4
ifort = -O3 -traceback
#-qopenmp 
# Define F90 flags #
# -Mipa=fast is suppressed #
F90FLAGS= -cpp $(PGI) $(ifort) $(MULTICORE) $(COMPILE) $(DEBUG) $(DIVB) $(FIXPOLE) $(MATHPACK)
LDFLAGS= -cpp $(PGI) $(ifort) $(MULTICORE) $(COMPILE) $(DEBUG) $(DIVB) $(FIXPOLE) $(MATHPACK)
