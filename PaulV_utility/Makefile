# Makefile for PaulV utilitites
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

THIS = PaulV_utility

SRC_NAMES = Type_Kinds File_Utility Error_Handler Fundamental_Constants Endian_Utility Binary_File_Utility \
            String_Utility Compare_Float_Numbers Date_Utility List_File_Utility
						
SRC := $(addsuffix .f90, $(SRC_NAMES))
OBJ := $(addsuffix .o, $(SRC_NAMES))
MOD := $(addsuffix .mod, $(SRC_NAMES))


# Object recipies

all : $(THIS)

$(OBJ) : $(SRC)
	$(FC) -c $(FOPT) $(addsuffix .f90, $(basename $@))
	
$(MOD) : $(OBJ)
	cp -f $@ ../mods/$(PLATFORM)/


# Final recipies

$(THIS) : $(MOD)
	$(AR) cr lib$@.a $(OBJ)
	cp -f lib$@.a ../libs/$(PLATFORM)/

clean : 
	rm -rf *.o *.mod *.a

	
