#=====================================================================
# tovSolve, C version
#
#   make          build tov_c.x (radius) and tov_h_c.x (enthalpy)
#   make clean    remove objects and executables
#
# Run from the repository root, where the EOS tables are, e.g.
#   c/tov_h_c.x eos_SLY4.in 0.3 1.2 4
#=====================================================================
CC     = gcc
CFLAGS = -O2 -Wall -Wextra -std=c99
LDLIBS = -lm

PROGS = tov_c.x tov_h_c.x

.PHONY : all clean

all : $(PROGS)

tov_c.x : tov.o tov_main.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

tov_h_c.x : tov.o tov_h.o tov_main_h.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

# Header dependencies (tov_main_h.c includes tov_main.c)
tov.o tov_h.o             : tov.h tov_internal.h
tov_main.o                : tov.h
tov_main_h.o              : tov.h tov_main.c

clean :
	rm -f *.o $(PROGS)
