#ifndef TOV_H
#define TOV_H

#include <stdio.h>

#define TOV_FM3CM3 1.0e39

/* EOS table, RNS format. Stored as log10 of the file columns. */
typedef struct {
    int np;        /* Number of rows */
    double *eray;  /* Energy density/c^2 (g/cm^3) */
    double *pray;  /* Pressure (dyn/cm^2) */
    double *hray;  /* Enthalpy (cm^2/s^2), not used by the solvers */
    double *xnray; /* Baryon number density (1/cm^3) */
} TovEos;

typedef struct {
    double mass;           /* M (Msun) */
    double radius;         /* R (km) */
    double k2;             /* Tidal Love number */
    double lambda;         /* Tidal deformability (1e36 g cm^2 s^2) */
    double moment_inertia; /* I (1e45 g cm^2) */
    double beta;           /* Compactness GM/(Rc^2) */
    double rhoc;           /* Central number density (fm^-3) */
} TovResult;

int tov_load_eos(const char *eos_file, TovEos *eos);
void tov_free_eos(TovEos *eos);

/* Radius as the independent variable (libtov.f90: SOLVE_TOV) */
int solve_tov(const TovEos *eos, double rhoc, TovResult *result);

/* Enthalpy formalism (libtov_h.f90: SOLVE_TOV_H) */
int solve_tov_h(const TovEos *eos, double rhoc, TovResult *result);

void tov_print_header(FILE *stream);
void tov_print_result(FILE *stream, const TovResult *result);

#endif
