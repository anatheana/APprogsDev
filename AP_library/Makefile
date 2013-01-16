# Makefile for AP_library
# Antti Penttilä, 2012

# Users can change these, but they should be already defined if using root makefile

ifndef $(ROOT_MAKE)

  PLATFORM = WIN64
  FC = gfortran
  CC = gcc
  AR = ar
  FOPT = -O2

endif


# Users should not change these

THIS = AP_library

SRC_NAMES = ErrorCapsule hashtable AP_quicksort AP_utils InputParser
						
SRC := $(addsuffix .f95, $(SRC_NAMES))
OBJ := $(addsuffix .o, $(SRC_NAMES))
MOD := $(addsuffix .mod, $(SRC_NAMES))

INC_DIR := ../mods/$(PLATFORM)/

# Object recipies

all : $(THIS)

$(OBJ) : $(SRC)
	$(FC) -c $(FOPT) -I$(INC_DIR) $(addsuffix .f95, $(basename $@))
	
$(MOD) : $(OBJ)
	cp -f $@ ../mods/$(PLATFORM)/


# Final recipies

$(THIS) : $(MOD)
	$(AR) cr lib$@.a $(OBJ)
	cp -f lib$@.a ../libs/$(PLATFORM)/

clean : 
	rm -rf *.o *.mod *.a

	
