#############################################################
#
# Alexander Soudackov, The Pennsylvania State University
# August 5, 2005
#
#############################################################

.SUFFIXES: .o .mod .f90 .f

FC      = gfortran
FCFLAGS = -O2
PROGRAM = fgh_bspline
LDFLAGS =
EXEDIR  = ../bin
EXE     = $(EXEDIR)/$(PROGRAM).bin

# rules

%.o : %.f90
	$(FC) -c $(FCFLAGS) -o $*.o $*.f90

%.o : %.f
	$(FC) -c $(FCFLAGS) -o $*.o $*.f

%.mod : %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi
	
SRCOBJ = fgh_bspline.o schroedinger.o rs.o
SRCSRC = fgh_bspline.f90 schroedinger.f90 rs.f
MODOBJ = bspline_mod.o
MODMOD = bspline_mod.mod
MODSRC = bspline_mod.f90

build:  $(MODMOD)
	@echo "*** building object files ***"
	$(MAKE) FC=$(FC) $(SRCOBJ) $(MODOBJ)

all:    build
	@echo "*** linking $(PROGRAM) ***"
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $(EXE) $(SRCOBJ) $(MODOBJ)
	@echo "---------------------------------------------------------------------------------------"
	@size $(EXE)
	@echo "---------------------------------------------------------------------------------------"

clean:
	rm -f *.o *.mod

#############################
# source files dependencies
#############################

fgh_bspline.o : fgh_bspline.f90 bspline_mod.o

#############################
# module files dependencies
#############################
