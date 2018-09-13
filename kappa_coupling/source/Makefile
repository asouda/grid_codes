#############################################################
#
# Alexander Soudackov, The Pennsylvania State University
# August 5, 2005
#
#############################################################

.SUFFIXES: .o .mod .f90 .f

FC      = gfortran
FCFLAGS = -O2
PROGRAM = kappa_coupling
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
	
SRCOBJ = kappa_coupling.o schroedinger.o rs.o dgamma.o
SRCSRC = kappa_coupling.f90 schroedinger.f90 rs.f dgamma.f
MODOBJ = module_timers.o module_bspline.o
MODMOD = timers.mod bspline.mod numeric.mod
MODSRC = module_timers.f90 module_bspline.f90

build:
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

kappa_coupling.o : kappa_coupling.f90 module_timers.o module_bspline.o rs.o schroedinger.o dgamma.o

#############################
# module files dependencies
#############################
