/*
 * TOV solver in the enthalpy formalism
 *
 * L. Lindblom, ApJ 398, 569 (1992)
 *
 * C version of libtov_h.f90. The independent variable is the
 * pseudo-enthalpy
 *    h = int_0^P dP' / ( e(P') + P' )    ( c = G = 1 )
 * measured from the surface, so the star is integrated from h = h_c at
 * the center to exactly h = 0 at the surface. The surface is the same as
 * in solve_tov: P = p_term (second row of the EOS table).
 *
 * P is carried as a dependent variable ( dP/dh = e + P ), and e(P),
 * cs2(P) come from tov_eos_interp, so this solves the same EOS as
 * solve_tov. The enthalpy column of the EOS file is not used.
 */

#include "tov_internal.h"

#include <math.h>
#include <stdio.h>

enum {
    TOV_H_MAX_STEPS = 100000
};

/*=====================================================================*/
/* Enthalpy difference between pressures p1 < p2                       */
/*    delta_h = int_p1^p2 dP / ( e(P) c^2 + P )    (dimensionless)     */
/* with e(P) from tov_eos_interp. The interpolant is smooth between    */
/* table points, so 5-point Gauss-Legendre in log10(P) on each table   */
/* interval is accurate to round-off                                   */
/*=====================================================================*/
static double delta_h(const TovEos *eos, double p1, double p2)
{
    const double ln10 = 2.302585092994045684;
    const double tg[5] = {-0.9061798459386640, -0.5384693101056831, 0.0,
                          0.5384693101056831, 0.9061798459386640};
    const double wg[5] = {0.2369268850561891, 0.4786286704993665,
                          0.5688888888888889, 0.4786286704993665,
                          0.2369268850561891};
    double xa = log10(p1);
    double xb = log10(p2);
    double s = 0.0;
    double a = xa;

    /* Split [xa, xb] at the table points (and past the end of the table) */
    for (int k = 0; k <= eos->np; k++) {
        double b;
        if (k < eos->np) {
            if (eos->pray[k] <= a) {
                continue;
            }
            b = fmin(eos->pray[k], xb);
        } else {
            b = xb;
        }

        for (int j = 0; j < 5; j++) {
            double x = 0.5 * (a + b) + 0.5 * (b - a) * tg[j];
            double p = pow(10.0, x);
            double eden;
            double cs2;
            tov_eos_interp(eos, p, &eden, &cs2);
            s = s + 0.5 * (b - a) * wg[j] * p * ln10 /
                        (eden * TOV_C * TOV_C + p);
        }

        a = b;
        if (a >= xb) {
            break;
        }
    }

    return s;
}

/*=====================================================================*/
/* TOV, y and f equations with the enthalpy as the independent         */
/* variable ( c = G = 1, lengths in m )                                */
/* y[0]: radius, y[1]: gravitational mass, y[2]: pressure,             */
/* y[3]: y = r H'/H (l = 2), y[4]: f = d ln(omega) / d ln(r)           */
/* The equations do not depend on the enthalpy explicitly              */
/*=====================================================================*/
static void tov_h_derivs(double hent, const double y[], double dydh[],
                         const void *ctx)
{
    const TovEos *eos = ctx;
    const double pi = TOV_PI;
    double r = y[0];
    double m = y[1];
    double p = y[2];
    double yy = y[3];
    double ff = y[4];
    double eden;
    double ed;
    double cs2;
    double drdh;
    double r2, r3, r4, m2;
    double f;
    double q;
    double l;
    double s;

    (void)hent;

    tov_eos_interp(eos, p / TOV_PG, &eden, &cs2);
    ed = eden * TOV_C * TOV_C * TOV_PG;

    r2 = r * r;
    r3 = r2 * r;
    r4 = r2 * r2;
    m2 = m * m;
    l = 1.0 / (1.0 - 2.0 * m / r);

    /* Radius, mass and pressure */
    drdh = -r * (r - 2.0 * m) / (m + 4.0 * pi * r3 * p);
    dydh[0] = drdh;
    dydh[1] = (4.0 * pi) * r2 * ed * drdh;
    dydh[2] = ed + p;

    /* y: dy/dh = (dy/dr)(dr/dh) */
    f = (1.0 - (4.0 * pi * r2) * (ed - p)) * l;
    s = 1.0 + 4.0 * pi * r3 * p / m;
    q = (4.0 * pi) * (5.0 * ed + 9.0 * p + ((ed + p) / cs2)) * l -
        (6.0 / r2) * l - (4.0 * m2 / r4) * (s * s) * (l * l);
    dydh[3] = (-(yy * yy / r) - (yy * f / r) - (r * q)) * drdh;

    /* f: df/dh = (df/dr)(dr/dh) */
    dydh[4] = (-(ff / r) * (ff + 3.0) +
               (4.0 + ff) * (4.0 * pi * r2) * (ed + p) / (r - 2.0 * m)) *
              drdh;
}

/*=====================================================================*/
/* Integrate y[] from x1 to x2 with tov_rkdp (x2 < x1 allowed)         */
/* h1, hmin: first and smallest allowed step size (magnitudes)         */
/*=====================================================================*/
static int hint(double y[], int nvar, double x1, double x2, double eps,
                double h1, double hmin, TovDerivFunc derivs,
                const void *ctx)
{
    const double tiny = 1.0e-300; /* P at the surface is ~1e-31 (1/m^2) */
    double x = x1;
    double h = copysign(h1, x2 - x1);
    double yscal[TOV_MAX_EQNS];
    double dydx[TOV_MAX_EQNS];

    if (nvar > TOV_MAX_EQNS) {
        fprintf(stderr, "too many equations in hint\n");
        return -1;
    }

    derivs(x, y, dydx, ctx);

    for (int nstp = 0; nstp < TOV_H_MAX_STEPS; nstp++) {
        double hdid;
        double hnext;

        for (int i = 0; i < nvar; i++) {
            yscal[i] = fabs(y[i]) + tiny;
        }

        /* Do not step past x2 */
        if ((x + h - x2) * (x + h - x1) > 0.0) {
            h = x2 - x;
        }

        if (tov_rkdp(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext,
                     derivs, ctx) != 0) {
            return -1;
        }

        if ((x - x2) * (x2 - x1) >= 0.0) {
            return 0;
        }

        h = hnext;
        if (fabs(h) < hmin) {
            fprintf(stderr,
                    "WARNING hint: step size below hmin before the surface "
                    "at x = %.15g\n",
                    x);
            return 0;
        }
    }

    fprintf(stderr,
            "WARNING hint: too many steps before the surface at x = %.15g\n",
            x);
    return 0;
}

/*=====================================================================*/
/* Compute spherical NS model {M, R, k2, lambda, I, beta, rhoc}        */
/* rhoc: central number density (1/cm^3)                              */
/*=====================================================================*/
int solve_tov_h(const TovEos *eos, double rhoc, TovResult *result)
{
    const double pi = TOV_PI;
    enum { nvar = 5 }; /* r, m, P, y, f */
    double p_c;   /* Central pressure (dyn/cm^2) */
    double edenc; /* Central energy density (g/cm^3) */
    double cs2c;  /* Central speed of sound squared */
    double e_g;   /* Central energy density (1/m^2) */
    double p_g;   /* Central pressure (1/m^2) */
    double dedh;  /* de/dh at the center (1/m^2) */
    double h_c;   /* Central enthalpy */
    double dh;    /* Offset of the starting point */
    double h0;    /* Starting enthalpy */
    double r0;    /* Starting radius (m) */
    double m0;    /* Starting mass (m) */
    double y[nvar];
    double p_term;
    const double eps = 1.0e-8;

    /* Central pressure, energy density and speed of sound */
    p_c = tov_eosinv(rhoc, eos->xnray, eos->pray, eos->np);
    tov_eos_interp(eos, p_c, &edenc, &cs2c);

    /* Termination pressure (surface, h = 0) and central enthalpy */
    p_term = tov_p_term(eos);
    h_c = delta_h(eos, p_term, p_c);

    /* Central values in geometric units ( c = G = 1, lengths in m ) */
    e_g = edenc * TOV_C * TOV_C * TOV_PG;
    p_g = p_c * TOV_PG;
    dedh = (e_g + p_g) / cs2c;

    /* Initial conditions at h0 = h_c - dh, from the series expansion  */
    /* around the center (Lindblom 1992, Eqs. (7) and (8))             */
    dh = 1.0e-7 * h_c;
    h0 = h_c - dh;
    r0 = sqrt(3.0 * dh / (2.0 * pi * (e_g + 3.0 * p_g))) *
         (1.0 - 0.25 * (e_g - 3.0 * p_g - 0.6 * dedh) * dh /
                    (e_g + 3.0 * p_g));
    m0 = (4.0 / 3.0) * pi * e_g * (r0 * r0 * r0) *
         (1.0 - 0.6 * dedh * dh / e_g);

    y[0] = r0;                                      /* Radius (m) */
    y[1] = m0;                                      /* Gravitational mass (m) */
    y[2] = p_g - (e_g + p_g) * dh;                  /* Pressure (1/m^2) */
    y[3] = 2.0;                                     /* y (for k2) */
    y[4] = (16.0 * pi / 5.0) * (e_g + p_g) * (r0 * r0); /* f (for I) */

    /* Integrate from the center (h0) to the surface (h = 0) */
    if (hint(y, nvar, h0, 0.0, eps, dh, 1.0e-14 * h_c, tov_h_derivs, eos) !=
        0) {
        return -1;
    }

    tov_star_output(y[1] / TOV_MG, y[0] * 1.0e2, y[3], y[4], rhoc, result);
    return 0;
}
