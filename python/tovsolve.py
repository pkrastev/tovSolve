#!/usr/bin/env python3
"""
tovsolve: neutron-star structure in the enthalpy formalism.

Python version of libtov_h.f90 / c/tov_h.c. For each central density it
integrates the TOV equations with the pseudo-enthalpy as the independent
variable (L. Lindblom, ApJ 398, 569 (1992)), together with the equations for
the tidal Love number k2 and the moment of inertia I, and returns
{M, R, k2, lambda, I, beta, rhoc}.

Every expression keeps the evaluation order of the Fortran and C versions,
so the results agree with them to the last bit. The solver is compiled with
Numba when it is installed (about C speed); otherwise it runs as plain
Python, about 100 times slower. Set TOVSOLVE_NO_NUMBA=1 to force plain Python.

Command line (same arguments and output as c/tov_h_c.x):
    python tovsolve.py [eos_file] [rho_start] [rho_end] [nsteps]

Library:
    import tovsolve
    eos = tovsolve.Eos("eos_SLY4.in")
    star = tovsolve.solve_tov_h(eos, 0.5)          # n_c in fm^-3
    seq = tovsolve.sequence(eos, 0.1, 1.5, 281)     # dict of numpy arrays
"""

import argparse
import math
import os
import sys
from collections import namedtuple

import numpy as np

if os.environ.get("TOVSOLVE_NO_NUMBA"):
    NUMBA = False
else:
    try:
        from numba import njit
        NUMBA = True
    except ImportError:
        NUMBA = False

if not NUMBA:
    def njit(*args, **kwargs):
        """No-op stand-in for numba.njit."""
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return lambda f: f

__all__ = ["Eos", "Star", "solve_tov_h", "sequence", "format_header",
           "format_row", "NUMBA"]

# Constants, as in MODULE CONSTANTS (libtov.f90)
PI = 3.1415926535897932384
C = 2.99792458e10
G = 6.67259e-8
M_SUN = 1.9892e33
FM3CM3 = 1.0e39
PG = (G / (C * C * C * C)) * 1.0e4  # e, P in 1/m^2
MG = (G / (C * C)) * 1.0e-2         # mass in m

EPS = 1.0e-8          # Relative accuracy of the integration
MAX_STEPS = 100000
TINY = 1.0e-300       # P at the surface is ~1e-31 (1/m^2)

Star = namedtuple("Star", "mass radius k2 lam I beta rhoc")
Star.__doc__ = """One neutron star: mass (Msun), radius (km), Love number k2,
tidal deformability lam (1e36 g cm^2 s^2), moment of inertia I (1e45 g cm^2),
compactness beta, central baryon density rhoc (fm^-3)."""


# ---------------------------------------------------------------------------
# EOS table
# ---------------------------------------------------------------------------
class Eos:
    """EOS table in RNS format, stored as log10 of the file columns.

    First line: number of rows; then per row energy density/c^2 (g/cm^3),
    pressure (dyn/cm^2), enthalpy (cm^2/s^2), baryon number density (cm^-3).
    The enthalpy column is read but not used by the solver.
    """

    def __init__(self, path):
        self.path = path
        try:
            with open(path) as f:
                tokens = f.read().split()
        except OSError as err:
            raise OSError(f"cannot open EOS file {path}: {err.strerror}") from err
        try:
            n = int(tokens[0])
        except (IndexError, ValueError):
            n = 0
        if n < 3:
            raise ValueError(f"bad number of rows in {path}")
        rows = []
        for i in range(n):
            try:
                row = [float(v) for v in tokens[1 + 4 * i:5 + 4 * i]]
            except ValueError:
                row = []
            if len(row) != 4 or min(row) <= 0.0:
                raise ValueError(f"cannot read row {i + 1} of {path}")
            rows.append(row)
        # math.log10 (the C library log10), not np.log10: NumPy's SIMD log10
        # can differ in the last bit, which would break agreement with C.
        a = np.array([[math.log10(v) for v in row] for row in rows])
        self.np = n
        self.eray = a[:, 0].copy()   # Energy density
        self.pray = a[:, 1].copy()   # Pressure
        self.hray = a[:, 2].copy()   # Enthalpy (not used)
        self.xnray = a[:, 3].copy()  # Baryon number density

    @property
    def n_max(self):
        """Largest baryon density in the table (fm^-3)."""
        return 10.0 ** self.xnray[-1] / FM3CM3

    def __repr__(self):
        name = os.path.basename(str(self.path))
        return f"Eos({name!r}, rows={self.np}, n_max={self.n_max:.3f} fm^-3)"


# ---------------------------------------------------------------------------
# Interpolation (as LAGINT / EOSINV / EOS_INTERP)
# ---------------------------------------------------------------------------
@njit(cache=True)
def _eosinv(xx, xi, yi):
    """3-point Lagrange interpolation in log10 space; returns 10**y."""
    ni = xi.shape[0]
    if xx <= 0.0:
        return 10.0 ** yi[0]
    x = math.log10(xx)
    if x <= xi[0]:
        return 10.0 ** yi[0]
    if x >= xi[ni - 1]:
        return 10.0 ** yi[ni - 1]
    n = 3
    i = 0
    j = ni - 1
    while j > i + 1:
        k = (i + j) // 2
        if x < xi[k]:
            j = k
        else:
            i = k
    i = i + 1 - n // 2
    if i < 0:
        i = 0
    if i + n > ni:
        i = ni - n
    y = 0.0
    for js in range(i, i + n):
        lam = 1.0
        for jl in range(i, i + n):
            if jl != js:
                lam = lam * (x - xi[jl]) / (xi[js] - xi[jl])
        y = y + yi[js] * lam
    return 10.0 ** y


@njit(cache=True)
def _eos_interp(p, pray, eray):
    """Energy density (g/cm^3) and cs2 = dP/dE (dimensionless) at pressure p."""
    np_ = pray.shape[0]
    x = math.log10(p) if p > 0.0 else pray[0]
    i = 0
    j = np_ - 1
    while j > i + 1:
        k = (i + j) // 2
        if x < pray[k]:
            j = k
        else:
            i = k
    if i + 3 > np_:
        i = np_ - 3
    x = min(max(x, pray[0]), pray[np_ - 1])
    x0 = pray[i]
    x1 = pray[i + 1]
    x2 = pray[i + 2]
    w0 = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2))
    w1 = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))
    w2 = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1))
    dw0 = ((x - x1) + (x - x2)) / ((x0 - x1) * (x0 - x2))
    dw1 = ((x - x0) + (x - x2)) / ((x1 - x0) * (x1 - x2))
    dw2 = ((x - x0) + (x - x1)) / ((x2 - x0) * (x2 - x1))
    eden = 10.0 ** (w0 * eray[i] + w1 * eray[i + 1] + w2 * eray[i + 2])
    dle = dw0 * eray[i] + dw1 * eray[i + 1] + dw2 * eray[i + 2]
    cs2 = (10.0 ** x / (eden * C * C)) / dle
    return eden, cs2


@njit(cache=True)
def _delta_h(p1, p2, pray, eray):
    """int_p1^p2 dP / (e(P) c^2 + P), 5-point Gauss-Legendre per table interval."""
    ln10 = 2.302585092994045684
    tg = (-0.9061798459386640, -0.5384693101056831, 0.0,
          0.5384693101056831, 0.9061798459386640)
    wg = (0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
          0.4786286704993665, 0.2369268850561891)
    np_ = pray.shape[0]
    xa = math.log10(p1)
    xb = math.log10(p2)
    s = 0.0
    a = xa
    for k in range(np_ + 1):
        if k < np_:
            if pray[k] <= a:
                continue
            b = min(pray[k], xb)
        else:
            b = xb
        for j in range(5):
            x = 0.5 * (a + b) + 0.5 * (b - a) * tg[j]
            p = 10.0 ** x
            eden, cs2 = _eos_interp(p, pray, eray)
            s = s + 0.5 * (b - a) * wg[j] * p * ln10 / (eden * C * C + p)
        a = b
        if a >= xb:
            break
    return s


# ---------------------------------------------------------------------------
# Equations and integrator
# ---------------------------------------------------------------------------
@njit(cache=True)
def _derivs(y, dydh, pray, eray):
    """TOV, y and f equations with the enthalpy as the independent variable
    (c = G = 1, lengths in m). y = [r, m, P, y, f]."""
    pi = PI
    r = y[0]
    m = y[1]
    p = y[2]
    yy = y[3]
    ff = y[4]
    eden, cs2 = _eos_interp(p / PG, pray, eray)
    ed = eden * C * C * PG
    r2 = r * r
    r3 = r2 * r
    r4 = r2 * r2
    m2 = m * m
    l = 1.0 / (1.0 - 2.0 * m / r)

    # Radius, mass and pressure
    drdh = -r * (r - 2.0 * m) / (m + 4.0 * pi * r3 * p)
    dydh[0] = drdh
    dydh[1] = (4.0 * pi) * r2 * ed * drdh
    dydh[2] = ed + p

    # y: dy/dh = (dy/dr)(dr/dh)
    f = (1.0 - (4.0 * pi * r2) * (ed - p)) * l
    s = 1.0 + 4.0 * pi * r3 * p / m
    q = ((4.0 * pi) * (5.0 * ed + 9.0 * p + ((ed + p) / cs2)) * l
         - (6.0 / r2) * l - (4.0 * m2 / r4) * (s * s) * (l * l))
    dydh[3] = (-(yy * yy / r) - (yy * f / r) - (r * q)) * drdh

    # f: df/dh = (df/dr)(dr/dh)
    dydh[4] = (-(ff / r) * (ff + 3.0)
               + (4.0 + ff) * (4.0 * pi * r2) * (ed + p) / (r - 2.0 * m)) * drdh


@njit(cache=True)
def _rkdp(y, dydx, x, htry, eps, yscal, work, pray, eray):
    """Adaptive Dormand-Prince 5(4) step (J.R. Dormand & P.J. Prince,
    J. Comp. Appl. Math. 6, 19 (1980)). Updates y and dydx in place (FSAL).
    Returns (x_new, hdid, hnext, status); status -1 = step size underflow."""
    c2 = 1.0 / 5.0
    c3 = 3.0 / 10.0
    c4 = 4.0 / 5.0
    c5 = 8.0 / 9.0
    a21 = 1.0 / 5.0
    a31 = 3.0 / 40.0
    a32 = 9.0 / 40.0
    a41 = 44.0 / 45.0
    a42 = -56.0 / 15.0
    a43 = 32.0 / 9.0
    a51 = 19372.0 / 6561.0
    a52 = -25360.0 / 2187.0
    a53 = 64448.0 / 6561.0
    a54 = -212.0 / 729.0
    a61 = 9017.0 / 3168.0
    a62 = -355.0 / 33.0
    a63 = 46732.0 / 5247.0
    a64 = 49.0 / 176.0
    a65 = -5103.0 / 18656.0
    b1 = 35.0 / 384.0          # 5th order weights (also the 7th stage row)
    b3 = 500.0 / 1113.0
    b4 = 125.0 / 192.0
    b5 = -2187.0 / 6784.0
    b6 = 11.0 / 84.0
    e1 = 71.0 / 57600.0        # Error weights: 5th order minus 4th order
    e3 = -71.0 / 16695.0
    e4 = 71.0 / 1920.0
    e5 = -17253.0 / 339200.0
    e6 = 22.0 / 525.0
    e7 = -1.0 / 40.0
    safety = 0.9
    pgrow = -0.2
    pshrnk = -0.25
    errcon = 1.89e-4           # (5/safety)**(1/pgrow)

    n = y.shape[0]
    k2 = work[0]
    k3 = work[1]
    k4 = work[2]
    k5 = work[3]
    k6 = work[4]
    k7 = work[5]
    ytemp = work[6]
    h = htry
    while True:
        for i in range(n):
            ytemp[i] = y[i] + h * a21 * dydx[i]
        _derivs(ytemp, k2, pray, eray)
        for i in range(n):
            ytemp[i] = y[i] + h * (a31 * dydx[i] + a32 * k2[i])
        _derivs(ytemp, k3, pray, eray)
        for i in range(n):
            ytemp[i] = y[i] + h * (a41 * dydx[i] + a42 * k2[i] + a43 * k3[i])
        _derivs(ytemp, k4, pray, eray)
        for i in range(n):
            ytemp[i] = y[i] + h * (a51 * dydx[i] + a52 * k2[i] + a53 * k3[i]
                                   + a54 * k4[i])
        _derivs(ytemp, k5, pray, eray)
        for i in range(n):
            ytemp[i] = y[i] + h * (a61 * dydx[i] + a62 * k2[i] + a63 * k3[i]
                                   + a64 * k4[i] + a65 * k5[i])
        xnew = x + h
        _derivs(ytemp, k6, pray, eray)
        for i in range(n):
            ytemp[i] = y[i] + h * (b1 * dydx[i] + b3 * k3[i] + b4 * k4[i]
                                   + b5 * k5[i] + b6 * k6[i])
        _derivs(ytemp, k7, pray, eray)

        errmax = 0.0
        for i in range(n):
            yerr = h * (e1 * dydx[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i]
                        + e6 * k6[i] + e7 * k7[i])
            errmax = max(errmax, abs(yerr / yscal[i]))
        errmax = errmax / eps
        if errmax <= 1.0:
            break
        # Reject: shrink the step, but by no more than a factor of 10
        h = math.copysign(max(abs(safety * h * errmax ** pshrnk),
                              0.1 * abs(h)), h)
        if x + h == x:
            return x, 0.0, 0.0, -1

    # Accept: grow the step, but by no more than a factor of 5
    if errmax > errcon:
        hnext = safety * h * errmax ** pgrow
    else:
        hnext = 5.0 * h
    for i in range(n):
        y[i] = ytemp[i]
        dydx[i] = k7[i]
    return xnew, h, hnext, 0


@njit(cache=True)
def _hint(y, x1, x2, eps, h1, hmin, pray, eray):
    """Integrate y from x1 to x2. Returns (status, x): 0 = reached x2,
    1 = step below hmin, 2 = too many steps, -1 = step size underflow."""
    n = y.shape[0]
    dydx = np.empty(n)
    yscal = np.empty(n)
    work = np.empty((7, n))
    x = x1
    h = math.copysign(h1, x2 - x1)
    _derivs(y, dydx, pray, eray)
    for nstp in range(MAX_STEPS):
        for i in range(n):
            yscal[i] = abs(y[i]) + TINY
        # Do not step past x2
        if (x + h - x2) * (x + h - x1) > 0.0:
            h = x2 - x
        x, hdid, hnext, status = _rkdp(y, dydx, x, h, eps, yscal, work,
                                       pray, eray)
        if status != 0:
            return status, x
        if (x - x2) * (x2 - x1) >= 0.0:
            return 0, x
        h = hnext
        if abs(h) < hmin:
            return 1, x
    return 2, x


@njit(cache=True)
def _solve_k2(yr, beta):
    """Love number k2 from y(R) and the compactness."""
    b2 = beta * beta
    b5 = (b2 * beta) * b2  # beta**5, as gfortran evaluates it
    t = 1.0 - 2.0 * beta
    num = (8.0 / 5.0) * b5 * (t * t) * (2.0 - yr + 2.0 * beta * (yr - 1.0))
    den = (2.0 * beta * (6.0 - 3.0 * yr + 3.0 * beta * (5.0 * yr - 8.0))
           + 4.0 * beta * beta * beta
           * (13.0 - 11.0 * yr + beta * (3.0 * yr - 2.0)
              + 2.0 * beta * beta * (1.0 + yr))
           + 3.0 * (t * t) * (2.0 - yr + 2.0 * beta * (yr - 1.0))
           * math.log(1.0 - 2.0 * beta))
    return num / den


@njit(cache=True)
def _solve(rhoc, pray, eray, xnray):
    """One star in the enthalpy formalism; rhoc in 1/cm^3.
    Returns (mass g, radius cm, yR, fR, status, x_end)."""
    pi = PI
    p_c = _eosinv(rhoc, xnray, pray)
    edenc, cs2c = _eos_interp(p_c, pray, eray)

    # Termination pressure (surface, h = 0) and central enthalpy
    eden_term = 10.0 ** eray[1]
    p_term = _eosinv(eden_term, eray, pray)
    h_c = _delta_h(p_term, p_c, pray, eray)

    # Central values in geometric units (c = G = 1, lengths in m)
    e_g = edenc * C * C * PG
    p_g = p_c * PG
    dedh = (e_g + p_g) / cs2c

    # Series expansion around the center (Lindblom 1992, Eqs. (7), (8))
    dh = 1.0e-7 * h_c
    h0 = h_c - dh
    r0 = (math.sqrt(3.0 * dh / (2.0 * pi * (e_g + 3.0 * p_g)))
          * (1.0 - 0.25 * (e_g - 3.0 * p_g - 0.6 * dedh) * dh
             / (e_g + 3.0 * p_g)))
    m0 = (4.0 / 3.0) * pi * e_g * (r0 * r0 * r0) * (1.0 - 0.6 * dedh * dh / e_g)

    y = np.empty(5)
    y[0] = r0                                        # Radius (m)
    y[1] = m0                                        # Gravitational mass (m)
    y[2] = p_g - (e_g + p_g) * dh                    # Pressure (1/m^2)
    y[3] = 2.0                                       # y (for k2)
    y[4] = (16.0 * pi / 5.0) * (e_g + p_g) * (r0 * r0)  # f (for I)

    status, x = _hint(y, h0, 0.0, EPS, dh, 1.0e-14 * h_c, pray, eray)
    return y[1] / MG, y[0] * 1.0e2, y[3], y[4], status, x


@njit(cache=True)
def _star_output(mass, radius, yr, fr, rhoc):
    """k2, lambda, I, beta from the surface values (as STAR_OUTPUT)."""
    g = G
    c = C
    r2 = radius * radius
    r3 = r2 * radius
    r5 = (r2 * radius) * r2  # R**5, as gfortran evaluates it
    beta = (g * mass) / (radius * c * c)
    k2 = _solve_k2(yr, beta)
    lam = (2.0 / (3.0 * g)) * k2 * r5
    i_ns = ((r3 * fr) / (6.0 + 2.0 * fr)) * (c * c / g)
    return (mass / M_SUN, radius / 1.0e5, k2, lam / 1.0e36, i_ns / 1.0e45,
            beta, rhoc / FM3CM3)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
_WARN = {1: "step size below hmin before the surface",
         2: "too many steps before the surface"}


def _solve_cm3(eos, rhoc):
    """Star for a central number density rhoc in 1/cm^3."""
    mass, radius, yr, fr, status, x = _solve(rhoc, eos.pray, eos.eray, eos.xnray)
    if status < 0:
        raise RuntimeError("stepsize underflow in rkdp")
    if status > 0:
        print(f"WARNING hint: {_WARN[status]} at x = {x:.15g}", file=sys.stderr)
    return Star(*_star_output(mass, radius, yr, fr, rhoc))


def solve_tov_h(eos, n_c):
    """Neutron star with central baryon density n_c (fm^-3). Returns a Star."""
    return _solve_cm3(eos, n_c * FM3CM3)


def sequence(eos, rho_start=0.09, rho_end=1.5, nsteps=100):
    """Sequence of stars for n_c = rho_start ... rho_end (fm^-3), with the same
    density steps as the Fortran and C drivers. Returns a dict of numpy arrays
    (mass, radius, k2, lam, I, beta, rhoc) plus the dimensionless Lambda."""
    if nsteps < 2:
        raise ValueError("nsteps must be at least 2")
    rho_h = (rho_end - rho_start) / float(nsteps - 1)
    rho_tmp = rho_start
    stars = []
    for _ in range(nsteps):
        stars.append(_solve_cm3(eos, rho_tmp * FM3CM3))
        rho_tmp = rho_tmp + rho_h
    out = {k: np.array([getattr(s, k) for s in stars]) for k in Star._fields}
    out["Lambda"] = (2.0 / 3.0) * out["k2"] / out["beta"] ** 5
    return out


def format_header():
    """Column header, as the Fortran and C programs print it."""
    return ("         M           R            k2         lambda          I"
            "          beta         rhoc")


def format_row(star):
    """One output row in the Fortran format 7(2x,f11.6)."""
    fields = []
    for v in star:
        s = f"{v:11.6f}"
        fields.append("  " + ("*" * 11 if len(s) > 11 else s))
    return "".join(fields)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Neutron-star sequence in the enthalpy formalism "
                    "(same arguments and output as c/tov_h_c.x).")
    ap.add_argument("eos_file", nargs="?", default="eos_MDI_x0.0.in",
                    help="EOS table in RNS format (default: %(default)s)")
    ap.add_argument("rho_start", nargs="?", type=float, default=0.09,
                    help="first central baryon density, fm^-3 (default: %(default)s)")
    ap.add_argument("rho_end", nargs="?", type=float, default=1.5,
                    help="last central baryon density, fm^-3 (default: %(default)s)")
    ap.add_argument("nsteps", nargs="?", type=int, default=100,
                    help="number of stars (default: %(default)s)")
    args = ap.parse_args(argv)
    if args.nsteps < 2:
        ap.error("nsteps must be at least 2")
    try:
        eos = Eos(args.eos_file)
    except (OSError, ValueError) as err:
        print(f"ERROR load_eos: {err}", file=sys.stderr)
        return 1

    print(format_header())
    rho_h = (args.rho_end - args.rho_start) / float(args.nsteps - 1)
    rho_tmp = args.rho_start
    for _ in range(args.nsteps):
        print(format_row(_solve_cm3(eos, rho_tmp * FM3CM3)))
        rho_tmp = rho_tmp + rho_h
    return 0


if __name__ == "__main__":
    sys.exit(main())
