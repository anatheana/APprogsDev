# Makefile for PackingProperties
# Antti Penttilä, 2013

# Users can change these, but they should be already defined if using root makefile

ifndef $(ROOT_MAKE)

  PLATFORM = WIN64
  FC = gfortran
  CC = gcc
  AR = ar
  FOPT = -O2

endif


# Users should not change these

THIS = packingproperties

SRC_NAMES = PackingProperties
						
SRC := $(addsuffix .f95, $(SRC_NAMES))
OBJ := $(addsuffix .o, $(SRC_NAMES))
PROG := $(addsuffix .exe, $(SRC_NAMES))

INC_DIR := ../mods/$(PLATFORM)/
LIB_DIR := ../libs/$(PLATFORM)/
BIN_DIR := ../bin/$(PLATFROM)/
LOCAL_LIBS := -lPaulV_utility -lAP_library


# Object recipies

all : $(THIS)

$(OBJ) : $(SRC)
	$(FC) -c -I$(INC_DIR) $(FOPT) $(addsuffix .f95, $(basename $@))
	
$(PROG) : $(OBJ)
	$(FC) -o $@ -I$(INC_DIR) -L$(LIB_DIR) $(FOPT) $(OBJ) $(LOCAL_LIBS) 
	cp -f $@ ../bin/$(PLATFORM)/


# Final recipies

$(THIS) : $(PROG)
	@cp -f $(PROG) ../bin/$(PLATFORM)/

clean : 
	rm -rf *.o *.exe

