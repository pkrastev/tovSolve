/*
 * TOV solver with the radius as the independent variable.
 * C version of libtov.f90: same equations, interpolation, stepper and
 * evaluation order, so it gives the same results as the Fortran code.
 */

#include "tov_internal.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    TOV_XDIM = 1000,     /* Number of rows for storage array (as XDIM) */
    TOV_MAX_STEPS = 20000
};

/*=====================================================================*/
/* Load EOS file                                                       */
/* Format: first line is the number of rows, then one row per point    */
/*         {Energy_Density, P, H, Number_Density} in CGS units (RNS)   */
/*=====================================================================*/
int tov_load_eos(const char *eos_file, TovEos *eos)
{
    FILE *fp;
    int np;

    memset(eos, 0, sizeof(*eos));
    fp = fopen(eos_file, "r");
    if (!fp) {
        fprintf(stderr, "ERROR load_eos: cannot open EOS file %s: %s\n",
                eos_file, strerror(errno));
        return -1;
    }

    if (fscanf(fp, "%d", &np) != 1 || np < 3) {
        fprintf(stderr, "ERROR load_eos: bad number of rows in %s\n",
                eos_file);
        fclose(fp);
        return -1;
    }

    eos->np = np;
    eos->eray = calloc((size_t)np, sizeof(*eos->eray));
    eos->pray = calloc((size_t)np, sizeof(*eos->pray));
    eos->hray = calloc((size_t)np, sizeof(*eos->hray));
    eos->xnray = calloc((size_t)np, sizeof(*eos->xnray));
    if (!eos->eray || !eos->pray || !eos->hray || !eos->xnray) {
        fprintf(stderr, "ERROR load_eos: cannot allocate table for %s\n",
                eos_file);
        fclose(fp);
        tov_free_eos(eos);
        return -1;
    }

    for (int i = 0; i < np; i++) {
        double eden;
        double p;
        double h0;
        double n0;
        if (fscanf(fp, "%lf %lf %lf %lf", &eden, &p, &h0, &n0) != 4 ||
            eden <= 0.0 || p <= 0.0 || h0 <= 0.0 || n0 <= 0.0) {
            fprintf(stderr, "ERROR load_eos: cannot read row %d of %s\n",
                    i + 1, eos_file);
            fclose(fp);
            tov_free_eos(eos);
            return -1;
        }
        eos->eray[i] = log10(eden);
        eos->pray[i] = log10(p);
        eos->hray[i] = log10(h0);
        eos->xnray[i] = log10(n0);
    }

    fclose(fp);
    return 0;
}

void tov_free_eos(TovEos *eos)
{
    free(eos->eray);
    free(eos->pray);
    free(eos->hray);
    free(eos->xnray);
    memset(eos, 0, sizeof(*eos));
}

/*=====================================================================*/
/* Energy density and speed of sound at pressure p                     */
/* eden -- energy density (g/cm^3), as EOSINV(p, pray, eray)           */
/* cs2  -- dP/dE with E = eden*c^2 (dimensionless)                     */
/* Same 3-point Lagrange interpolation in log10 space as LAGINT.       */
/* cs2 = (P/E) / (dlogE/dlogP), from the derivative of the interpolant */
/*=====================================================================*/
void tov_eos_interp(const TovEos *eos, double p, double *eden, double *cs2)
{
    const double *pray = eos->pray;
    const double *eray = eos->eray;
    int np = eos->np;
    int i;
    int j;
    double x;
    double x0;
    double x1;
    double x2;
    double w0, w1, w2;    /* Lagrange weights */
    double dw0, dw1, dw2; /* Derivatives of Lagrange weights */
    double dle;

    x = (p > 0.0) ? log10(p) : pray[0];

    /* Bisection: pray[i] <= x < pray[i+1] */
    i = 0;
    j = np - 1;
    while (j > i + 1) {
        int k = (i + j) / 2;
        if (x < pray[k]) {
            j = k;
        } else {
            i = k;
        }
    }
    if (i + 3 > np) {
        i = np - 3;
    }

    /* Outside the table take a boundary value, as in LAGINT */
    x = fmin(fmax(x, pray[0]), pray[np - 1]);

    x0 = pray[i];
    x1 = pray[i + 1];
    x2 = pray[i + 2];

    w0 = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2));
    w1 = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2));
    w2 = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
    dw0 = ((x - x1) + (x - x2)) / ((x0 - x1) * (x0 - x2));
    dw1 = ((x - x0) + (x - x2)) / ((x1 - x0) * (x1 - x2));
    dw2 = ((x - x0) + (x - x1)) / ((x2 - x0) * (x2 - x1));

    *eden = pow(10.0, w0 * eray[i] + w1 * eray[i + 1] + w2 * eray[i + 2]);

    dle = dw0 * eray[i] + dw1 * eray[i + 1] + dw2 * eray[i + 2];
    *cs2 = (pow(10.0, x) / (*eden * TOV_C * TOV_C)) / dle;
}

/*=====================================================================*/
/* Lagrange interpolation with n points (order n-1)                    */
/*=====================================================================*/
static double lagint(double xx, const double xi[], const double yi[], int ni,
                     int n)
{
    int i;
    int j;
    double y = 0.0;

    if (n > ni) {
        n = ni;
    }
    if (xx <= xi[0]) {
        return yi[0];
    }
    if (xx >= xi[ni - 1]) {
        return yi[ni - 1];
    }

    i = 0;
    j = ni - 1;
    while (j > i + 1) {
        int k = (i + j) / 2;
        if (xx < xi[k]) {
            j = k;
        } else {
            i = k;
        }
    }

    i = i + 1 - n / 2;
    if (i < 0) {
        i = 0;
    }
    if (i + n > ni) {
        i = ni - n;
    }

    for (int js = i; js < i + n; js++) {
        double lambda = 1.0;
        for (int jl = i; jl < i + n; jl++) {
            if (jl != js) {
                lambda = lambda * (xx - xi[jl]) / (xi[js] - xi[jl]);
            }
        }
        y = y + yi[js] * lambda;
    }

    return y;
}

/*=====================================================================*/
/* Interface to LAGINT: returns pressure, energy density, or number    */
/* density (tables in log10)                                           */
/*=====================================================================*/
double tov_eosinv(double xx, const double xi[], const double yi[], int ni)
{
    if (xx <= 0.0) {
        return pow(10.0, yi[0]);
    }
    return pow(10.0, lagint(log10(xx), xi, yi, ni, 3));
}

/* Termination pressure: second row of the EOS table */
double tov_p_term(const TovEos *eos)
{
    double eden_term = pow(10.0, eos->eray[1]);
    return tov_eosinv(eden_term, eos->eray, eos->pray, eos->np);
}

/*=====================================================================*/
/* TOV equations (cgs), together with the equations for y (tidal Love  */
/* number k2) and f (moment of inertia I)                              */
/* y[0]: gravitational mass, y[1]: pressure,                           */
/* y[2]: y = r H'/H (l = 2), y[3]: f = d ln(omega) / d ln(r)           */
/*=====================================================================*/
static void tov_derivs(double r, const double y[], double dydr[],
                       const void *ctx)
{
    const TovEos *eos = ctx;
    const double pi = TOV_PI;
    const double g = TOV_G;
    const double c2 = TOV_C * TOV_C;
    double grmass = y[0];
    double p = y[1];
    double yy = y[2];
    double ff = y[3];
    double eden;
    double cs2;
    double relcor;
    double r_g, m_g, ed_g, p_g;
    double r2, r3, r4, m2;
    double f;
    double q;
    double l;
    double s;

    tov_eos_interp(eos, p, &eden, &cs2);

    /* Gravitational mass */
    dydr[0] = (4.0 * pi) * r * r * eden;

    /* Pressure */
    relcor = (1.0 + p / (eden * c2)) *
             (1.0 + ((4.0 * pi) * p * r * r * r) / (grmass * c2)) /
             (1.0 - (2.0 * g * grmass) / (r * c2));
    dydr[1] = -(g * grmass / (r * r)) * eden * relcor;

    /* y and f in geometric units (c = G = 1, lengths in m) */
    r_g = r / 1.0e2;
    m_g = grmass * TOV_MG;
    ed_g = eden * c2 * TOV_PG;
    p_g = p * TOV_PG;
    r2 = r_g * r_g;
    r3 = r2 * r_g;
    r4 = r2 * r2;
    m2 = m_g * m_g;
    l = 1.0 / (1.0 - 2.0 * m_g / r_g);

    /* y */
    f = (1.0 - (4.0 * pi * r2) * (ed_g - p_g)) * l;
    s = 1.0 + 4.0 * pi * r3 * p_g / m_g;
    q = (4.0 * pi) * (5.0 * ed_g + 9.0 * p_g + ((ed_g + p_g) / cs2)) * l -
        (6.0 / r2) * l - (4.0 * m2 / r4) * (s * s) * (l * l);
    dydr[2] = (-(yy * yy / r_g) - (yy * f / r_g) - (r_g * q)) / 1.0e2;

    /* f */
    dydr[3] = (-(ff / r_g) * (ff + 3.0) +
               (4.0 + ff) * (4.0 * pi * r2) * (ed_g + p_g) /
                   (r_g - 2.0 * m_g)) /
              1.0e2;
}

/*=====================================================================*/
/* Evaluate Love number k2                                             */
/*=====================================================================*/
static double solve_k2(double yr, double beta)
{
    double b2 = beta * beta;
    double b5 = (b2 * beta) * b2; /* beta**5, as gfortran evaluates it */
    double t = 1.0 - 2.0 * beta;
    double num;
    double den;

    num = (8.0 / 5.0) * b5 * (t * t) * (2.0 - yr + 2.0 * beta * (yr - 1.0));
    den = 2.0 * beta * (6.0 - 3.0 * yr + 3.0 * beta * (5.0 * yr - 8.0)) +
          4.0 * beta * beta * beta *
              (13.0 - 11.0 * yr + beta * (3.0 * yr - 2.0) +
               2.0 * beta * beta * (1.0 + yr)) +
          3.0 * (t * t) * (2.0 - yr + 2.0 * beta * (yr - 1.0)) *
              log(1.0 - 2.0 * beta);
    return num / den;
}

/*=====================================================================*/
/* Compute k2, lambda, I, beta from the surface values                 */
/* mass (g), radius (cm), yr, fr at the surface, rhoc (1/cm^3)         */
/*=====================================================================*/
void tov_star_output(double mass, double radius, double yr, double fr,
                     double rhoc, TovResult *result)
{
    const double g = TOV_G;
    const double c = TOV_C;
    double r2 = radius * radius;
    double r3 = r2 * radius;
    double r5 = (r2 * radius) * r2; /* R**5, as gfortran evaluates it */
    double beta = (g * mass) / (radius * c * c);
    double k2 = solve_k2(yr, beta);
    double lambda = (2.0 / (3.0 * g)) * k2 * r5;
    double i_ns = ((r3 * fr) / (6.0 + 2.0 * fr)) * (c * c / g);

    result->mass = mass / TOV_M_SUN;
    result->radius = radius / 1.0e5;
    result->k2 = k2;
    result->lambda = lambda / 1.0e36;
    result->moment_inertia = i_ns / 1.0e45;
    result->beta = beta;
    result->rhoc = rhoc / TOV_FM3CM3;
}

/*=====================================================================*/
/* Adaptive Dormand-Prince 5(4) Runge-Kutta step                       */
/* Coefficients from J.R. Dormand & P.J. Prince,                       */
/* J. Comp. Appl. Math. 6, 19 (1980)                                   */
/* On return y, dydx and x are at x + hdid (FSAL: the last stage is    */
/* evaluated at the new point, so the driver does not call derivs)     */
/*=====================================================================*/
int tov_rkdp(double y[], double dydx[], int n, double *x, double htry,
             double eps, const double yscal[], double *hdid, double *hnext,
             TovDerivFunc derivs, const void *ctx)
{
    const double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0,
                 c5 = 8.0 / 9.0;
    const double a21 = 1.0 / 5.0;
    const double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    const double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    const double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0,
                 a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    const double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0,
                 a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0,
                 a65 = -5103.0 / 18656.0;
    /* 5th order weights (also the 7th stage row) */
    const double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0,
                 b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
    /* Error weights: 5th order minus 4th order */
    const double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0,
                 e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0,
                 e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;
    /* Step size control */
    const double safety = 0.9;
    const double pgrow = -0.2;
    const double pshrnk = -0.25;
    const double errcon = 1.89e-4; /* (5/safety)**(1/pgrow) */
    double k2[TOV_MAX_EQNS], k3[TOV_MAX_EQNS], k4[TOV_MAX_EQNS];
    double k5[TOV_MAX_EQNS], k6[TOV_MAX_EQNS], k7[TOV_MAX_EQNS];
    double ytemp[TOV_MAX_EQNS];
    double h = htry;
    double xnew;
    double errmax;

    if (n > TOV_MAX_EQNS) {
        fprintf(stderr, "too many equations in rkdp\n");
        return -1;
    }

    for (;;) {
        for (int i = 0; i < n; i++) ytemp[i] = y[i] + h * a21 * dydx[i];
        derivs(*x + c2 * h, ytemp, k2, ctx);
        for (int i = 0; i < n; i++)
            ytemp[i] = y[i] + h * (a31 * dydx[i] + a32 * k2[i]);
        derivs(*x + c3 * h, ytemp, k3, ctx);
        for (int i = 0; i < n; i++)
            ytemp[i] = y[i] + h * (a41 * dydx[i] + a42 * k2[i] + a43 * k3[i]);
        derivs(*x + c4 * h, ytemp, k4, ctx);
        for (int i = 0; i < n; i++)
            ytemp[i] = y[i] + h * (a51 * dydx[i] + a52 * k2[i] + a53 * k3[i] +
                                   a54 * k4[i]);
        derivs(*x + c5 * h, ytemp, k5, ctx);
        for (int i = 0; i < n; i++)
            ytemp[i] = y[i] + h * (a61 * dydx[i] + a62 * k2[i] + a63 * k3[i] +
                                   a64 * k4[i] + a65 * k5[i]);
        xnew = *x + h;
        derivs(xnew, ytemp, k6, ctx);
        for (int i = 0; i < n; i++)
            ytemp[i] = y[i] + h * (b1 * dydx[i] + b3 * k3[i] + b4 * k4[i] +
                                   b5 * k5[i] + b6 * k6[i]);
        derivs(xnew, ytemp, k7, ctx);

        errmax = 0.0;
        for (int i = 0; i < n; i++) {
            double yerr = h * (e1 * dydx[i] + e3 * k3[i] + e4 * k4[i] +
                               e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
            errmax = fmax(errmax, fabs(yerr / yscal[i]));
        }
        errmax = errmax / eps;

        if (errmax <= 1.0) {
            break;
        }

        /* Reject: shrink the step, but by no more than a factor of 10 */
        h = copysign(fmax(fabs(safety * h * pow(errmax, pshrnk)),
                          0.1 * fabs(h)),
                     h);
        if (*x + h == *x) {
            fprintf(stderr, "stepsize underflow in rkdp\n");
            return -1;
        }
    }

    /* Accept: grow the step, but by no more than a factor of 5 */
    if (errmax > errcon) {
        *hnext = safety * h * pow(errmax, pgrow);
    } else {
        *hnext = 5.0 * h;
    }

    *hdid = h;
    *x = xnew;
    for (int i = 0; i < n; i++) {
        y[i] = ytemp[i];
        dydx[i] = k7[i];
    }
    return 0;
}

/*=====================================================================*/
/* Integrate from x1 until the pressure y[1] drops to p_term (surface) */
/* Stores the steps as TOVINT does; on return xp[kount-1], yp[kount-1] */
/* hold the final point                                                */
/*=====================================================================*/
static int tovint(double ystart[], int nvar, double x1, double x2,
                  int kmax, int *kount, double xp[],
                  double yp[][TOV_MAX_EQNS], double eps, double h1,
                  double hmin, double p_term, TovDerivFunc derivs,
                  const void *ctx)
{
    const double tiny = 1.0e-15;
    double x = x1;
    double h = copysign(h1, x2 - x1);
    double y[TOV_MAX_EQNS];
    double dydx[TOV_MAX_EQNS];
    double yscal[TOV_MAX_EQNS];

    if (nvar > TOV_MAX_EQNS) {
        fprintf(stderr, "too many equations in tovint\n");
        return -1;
    }

    *kount = 0;
    for (int i = 0; i < nvar; i++) {
        y[i] = ystart[i];
    }

    /* Derivatives at the start; the stepper returns them at each new x */
    derivs(x, y, dydx, ctx);

    for (int nstp = 0; nstp < TOV_MAX_STEPS; nstp++) {
        double hdid;
        double hnext;

        for (int i = 0; i < nvar; i++) {
            yscal[i] = fabs(y[i]) + tiny;
        }

        /* Store intermediate results (every step, as dxsav = 0) */
        if (kmax > 0 && *kount < kmax - 1) {
            xp[*kount] = x;
            for (int i = 0; i < nvar; i++) {
                yp[*kount][i] = y[i];
            }
            (*kount)++;
        }

        /* If the step can overshoot the stop point then cut it */
        if ((x + h - x2) * (x + h - x1) > 0.0) {
            h = x2 - x;
        }

        if (tov_rkdp(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext,
                     derivs, ctx) != 0) {
            return -1;
        }

        /* Exit point: save the final step */
        if ((x - x2) * (x2 - x1) >= 0.0 || y[1] <= p_term) {
            for (int i = 0; i < nvar; i++) {
                ystart[i] = y[i];
            }
            if (kmax != 0) {
                xp[*kount] = x;
                for (int i = 0; i < nvar; i++) {
                    yp[*kount][i] = y[i];
                }
                (*kount)++;
            }
            return 0;
        }

        /* Set the step size for the next iteration */
        h = hnext;
        if (fabs(hnext) < hmin) {
            fprintf(stderr,
                    "WARNING tovint: step size below hmin before the "
                    "surface at x = %.15g\n",
                    x);
            return 0;
        }
    }

    fprintf(stderr,
            "WARNING tovint: too many steps before the surface at x = %.15g\n",
            x);
    return 0;
}

/*=====================================================================*/
/* Compute spherical NS model {M, R, k2, lambda, I, beta, rhoc}        */
/* rhoc: central number density (1/cm^3)                              */
/*=====================================================================*/
int solve_tov(const TovEos *eos, double rhoc, TovResult *result)
{
    const double pi = TOV_PI;
    const double g = TOV_G;
    enum { ydim = 4 }; /* M, P, y, f */
    double xp[TOV_XDIM];                 /* Storage for R */
    double yp[TOV_XDIM][TOV_MAX_EQNS];   /* Storage for M, P, y, f */
    double bc[ydim];
    double p_c;
    double edenc;
    double p_term;
    int kount;
    /* Integrate in radius (in cm) from h_start to h_stop */
    const double h_start = 1.0e2;
    const double h_try = h_start;
    const double h_min = 1.0e-6;
    const double h_stop = TOV_R_SUN;
    const double eps = 1.0e-8;

    /* Central pressure and energy density */
    p_c = tov_eosinv(rhoc, eos->xnray, eos->pray, eos->np);
    edenc = tov_eosinv(p_c, eos->pray, eos->eray, eos->np);

    /* Initial conditions */
    bc[0] = (4.0 / 3.0) * pi * (h_start * h_start * h_start) * edenc;
    bc[1] = p_c - 0.5 * (4.0 / 3.0 * pi) * g * (h_start * h_start) *
                      (edenc * edenc);
    bc[2] = 2.0;      /* y (for k2) */
    bc[3] = 0.000001; /* f (for I) */

    p_term = tov_p_term(eos);

    if (tovint(bc, ydim, h_start, h_stop, TOV_XDIM, &kount, xp, yp, eps,
               h_try, h_min, p_term, tov_derivs, eos) != 0) {
        return -1;
    }

    tov_star_output(yp[kount - 1][0], xp[kount - 1], yp[kount - 1][2],
                    yp[kount - 1][3], rhoc, result);
    return 0;
}

/*=====================================================================*/
/* Output in the Fortran format: header, then 7(2x,f11.6)              */
/*=====================================================================*/
void tov_print_header(FILE *stream)
{
    fprintf(stream,
            "         M           R            k2         lambda          I"
            "          beta         rhoc\n");
}

/* Like Fortran f11.6: 11 asterisks when the value does not fit */
static void print_f11_6(FILE *stream, double value)
{
    char buf[64];
    int len = snprintf(buf, sizeof(buf), "%11.6f", value);

    fprintf(stream, "  %s", (len > 11) ? "***********" : buf);
}

void tov_print_result(FILE *stream, const TovResult *result)
{
    print_f11_6(stream, result->mass);
    print_f11_6(stream, result->radius);
    print_f11_6(stream, result->k2);
    print_f11_6(stream, result->lambda);
    print_f11_6(stream, result->moment_inertia);
    print_f11_6(stream, result->beta);
    print_f11_6(stream, result->rhoc);
    fputc('\n', stream);
}
