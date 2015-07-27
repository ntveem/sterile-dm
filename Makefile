#-------------------------------------------------
# Makefile for 7keV 				 
#-------------------------------------------------

#-------------------------------------------------
# Environment variables -- Change if needed
#-------------------------------------------------
FC =gfortran

FFLAGS =-g -Wall -fbacktrace

# Library directory (if any)
LDIR = lib

# Subdirectory with object files
ODIR = obj

# compiler dependent flag to place module files in include directory,
# if possible 
MSFLAGS =-J$(ODIR)

# compiler dependent flag to load module files from include directory,
# if possible 
MLFLAGS =-I$(ODIR)

#-------------------------------------------------
# Dependent variables
#-------------------------------------------------
# list of all library files that will be included
LIBF = $(patsubst lib%.a, %, $(notdir $(wildcard $(LDIR)/lib*.a)))

# link any libraries needed 
ifeq ($(LIBF),)
	LFLAGS=-lm
else
	LFLAGS=-L$(LDIR) -lm $(addprefix -l,$(LIBF))
endif

#-------------------------------------------------
# Actions 
#-------------------------------------------------

include ./deps

clean: 
	rm -f sterile-nu $(ODIR)/*.o $(ODIR)/*.mod *.mod

#-------------------------------------------------
