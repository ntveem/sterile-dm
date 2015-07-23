#-------------------------------------------------
# Makefile for 7keV 				 
#-------------------------------------------------

#-------------------------------------------------
# Environment variables -- Change if needed
#-------------------------------------------------
FC =gfortran

FFLAGS =-O3 -ffast-math -funroll-loops -fopenmp

# Source directory
SRCDIR = src

# Subdirectory with module files, separate to enforce compilation order
MODDIR = src/modules

# Include directory
IDIR = include

# Library directory
LDIR = lib

# Subdirectory with object files
OBJDIR = obj

# compiler dependent flag to place module files in include directory,
# if possible 
MSFLAGS =-J$(IDIR)

# compiler dependent flag to load module files from include directory,
# if possible 
MLFLAGS =-I$(IDIR)

# purely to check that all the tables are present
DATADIR1 = data/tables
DATADIR2 = data/thermalpop

DATA1 = \
 $(DATADIR1)/ChiTable_alltemp.dat $(DATADIR1)/SMgstar.dat \
 $(DATADIR1)/dmudLe.dat $(DATADIR1)/dmudLmu.dat \
 $(DATADIR1)/dmudLtau.dat $(DATADIR1)/ParticleChiTable.dat \
 $(DATADIR1)/EFermi.dat $(DATADIR1)/rate_total.dat \

DATA2 = \
 $(DATADIR2)/thermmu.dat $(DATADIR2)/thermtau.dat \
 $(DATADIR2)/thermc.dat $(DATADIR2)/thermb.dat \

#-------------------------------------------------
# Dependent variables
#-------------------------------------------------
# list of all object files that will be created
OBJM = $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f, %.o, $(wildcard $(MODDIR)/*.f))))
OBJF = $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f, %.o, $(wildcard $(SRCDIR)/*.f))))

# list of all library files that will be included
LIBF = $(patsubst lib%.a, %, $(notdir $(wildcard $(LDIR)/lib*.a)))

# Include any other files from include directory 
IFLAGS =-I$(IDIR)

# link any libraries needed 
ifeq ($(LIBF),)
	LFLAGS=-lm
else
	LFLAGS=-L$(LDIR) -lm $(addprefix -l,$(LIBF))
endif

#-------------------------------------------------
# Actions 
#-------------------------------------------------
default: 7kev

all: 7kev

# link to get the main executable
7kev: $(OBJM) $(OBJF)
	$(FC) -o $@ $^ $(IFLAGS) $(LFLAGS)

# build module files first
$(OBJDIR)/%.o: $(MODDIR)/%.f $(DATA1) $(DATA2) 
	$(FC) $(FFLAGS) -c -o $@ $< $(MSFLAGS) $(IFLAGS) $(LFLAGS)

# then build other files
$(OBJDIR)/%.o: $(SRCDIR)/%.f
	$(FC) $(FFLAGS) -c -o $@ $< $(MLFLAGS) $(IFLAGS) $(LFLAGS)

clean: 
	rm -f 7kev $(OBJDIR)/*.o $(MODDIR)/*.mod *.mod

#-------------------------------------------------
