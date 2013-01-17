# Makefile for meshVolume
# Antti Penttilä, 2013

# Users can change these, but they should be already defined if using root makefile

ifndef $(ROOT_MAKE)

  PLATFORM = WIN64
  FC = gfortran
  CC = gcc
  AR = ar
  FOPT = -O2 -std=f2003 -ffree-form

endif


# Users should not change these

THIS = meshVolume

SRC_NAMES = meshVolume
						
SRC := $(addsuffix .f, $(SRC_NAMES))
OBJ := $(addsuffix .o, $(SRC_NAMES))
MOD := $(addsuffix .mod, $(SRC_NAMES))

INC_DIR := ../mods/$(PLATFORM)/

# Object recipies

all : $(THIS)

$(OBJ) : $(SRC)
	$(FC) -c $(FOPT) -I$(INC_DIR) $(addsuffix .f, $(basename $@))
	
$(MOD) : $(OBJ)
	cp -f $@ ../mods/$(PLATFORM)/


# Final recipies

$(THIS) : $(MOD)
	$(AR) cr lib$@.a $(OBJ)
	cp -f lib$@.a ../libs/$(PLATFORM)/

clean : 
	rm -rf *.o *.mod *.a
