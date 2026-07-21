#include header#
include ../../../make.inc

EXTRAINCS += $(HDF5INCS) -module .
EXTRAOBJECTS += $(HDF5LIBS)  -module .

SOLVER := $(CUMC3D_DIR)/solver/mhd
SOURCE := $(CUMC3D_DIR)/model/torus-cyn/src
CUSTOM := $(CUMC3D_DIR)/model/torus-cyn/custom-mhd

vpath %.f90 $(CUSTOM):$(SOLVER):$(SOURCE)
vpath %.h $(CUSTOM):$(SOLVER)

SOURCE1 = $(shell cat $(SOLVER)/source.list)
SOURCE2 = $(shell cat $(SOURCE)/source.list)

SRC = $(notdir $(SOURCE1)) $(notdir $(SOURCE2))
OBJ = $(patsubst %.f90,%.o,$(SRC))

ALL: $(OBJ)
	$(F90) $(LDFLAGS) -o ../CUMC3D $(OBJ) $(EXTRAOBJECTS)

# 7. THE WORKING PATTERN RULE
%.o: %.f90
	$(F90) $(F90FLAGS) $(EXTRAINCS) -c $< -o $@