#ifndef TOV_INTERNAL_H
#define TOV_INTERNAL_H

/* Shared by tov.c and tov_h.c. Not part of the public API (tov.h).     */
/* Expressions keep the evaluation order of the Fortran version, so the */
/* two give the same results to the last bit.                           */

#include "tov.h"

enum {
    TOV_MAX_EQNS = 20 /* Largest ODE system a driver or stepper accepts */
};

/* Constants, as in MODULE CONSTANTS (libtov.f90) */
#define TOV_PI 3.1415926535897932384
#define TOV_C 2.99792458e10
#define TOV_G 6.67259e-8
#define TOV_R_SUN 6.95997e10
#define TOV_M_SUN 1.9892e33
#define TOV_PG ((TOV_G / (TOV_C * TOV_C * TOV_C * TOV_C)) * 1.0e4) /* e, P in 1/m^2 */
#define TOV_MG ((TOV_G / (TOV_C * TOV_C)) * 1.0e-2)                /* mass in m */

typedef void (*TovDerivFunc)(double x, const double y[], double dydx[],
                             const void *ctx);

double tov_eosinv(double xx, const double xi[], const double yi[], int ni);
void tov_eos_interp(const TovEos *eos, double p, double *eden, double *cs2);
double tov_p_term(const TovEos *eos);
int tov_rkdp(double y[], double dydx[], int n, double *x, double htry,
             double eps, const double yscal[], double *hdid, double *hnext,
             TovDerivFunc derivs, const void *ctx);
void tov_star_output(double mass, double radius, double yr, double fr,
                     double rhoc, TovResult *result);

#endif
