#=====================================================================
# tovSolve, Fortran version
#
#   make          build tov.x (radius) and tov_h.x (enthalpy)
#   make clean    remove objects, modules and executables
#
# The C version has its own Makefile: make -C c
#=====================================================================
FC     = gfortran
FFLAGS = -O2

PROGS = tov.x tov_h.x

.PHONY : all clean

all : $(PROGS)

tov.x : libtov.o tov_main.o
	$(FC) $(FFLAGS) -o $@ $^

tov_h.x : libtov.o libtov_h.o tov_main_h.o
	$(FC) $(FFLAGS) -o $@ $^

%.o : %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

# These use the modules defined in libtov.f90
libtov_h.o tov_main.o tov_main_h.o : libtov.o

clean :
	rm -f *.o *.mod $(PROGS)
