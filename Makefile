#================================================
# Makefile
#================================================
F90COMPILER = gfortran
F90CFLAGS   = -c -O2
PRO         = tov

OBJECTS = libtov.o   \
          tov_main.o \
          ode.o

${PRO}.x : $(OBJECTS)
	$(F90COMPILER) -o ${PRO}.x $(OBJECTS)

%.o : %.f90
	$(F90COMPILER) $(F90CFLAGS) $(<F)

clean : 
	rm -rf *.o *.x *.mod
