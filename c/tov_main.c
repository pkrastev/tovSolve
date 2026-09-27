/*
 * Main program: sequence of neutron stars over a range of central
 * densities. tov_main_h.c builds the same driver with TOV_SOLVER set to
 * the enthalpy-formalism solver.
 *
 * Usage: prog [eos_file] [rho_start] [rho_end] [nsteps]
 */

#include "tov.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TOV_SOLVER
#define TOV_SOLVER solve_tov
#endif

static int parse_double(const char *text, const char *name, double *value)
{
    char *end = NULL;

    errno = 0;
    *value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0') {
        fprintf(stderr, "Invalid %s: '%s'\n", name, text);
        return -1;
    }
    return 0;
}

static int parse_int(const char *text, const char *name, int *value)
{
    char *end = NULL;
    long tmp;

    errno = 0;
    tmp = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || tmp <= 0) {
        fprintf(stderr, "Invalid %s: '%s'\n", name, text);
        return -1;
    }
    *value = (int)tmp;
    return 0;
}

int main(int argc, char **argv)
{
    const char *eos_file = "eos_MDI_x0.0.in"; /* EOS */
    double rho_start = 0.09;                  /* Starting density */
    double rho_end = 1.5;                     /* Ending density */
    int nsteps = 100;                         /* Number of density steps */
    double rho_h;
    double rho_tmp;
    TovEos eos;

    if (argc > 5) {
        fprintf(stderr,
                "Usage: %s [eos_file] [rho_start] [rho_end] [nsteps]\n",
                argv[0]);
        return EXIT_FAILURE;
    }
    if (argc > 1) {
        eos_file = argv[1];
    }
    if (argc > 2 && parse_double(argv[2], "rho_start", &rho_start) != 0) {
        return EXIT_FAILURE;
    }
    if (argc > 3 && parse_double(argv[3], "rho_end", &rho_end) != 0) {
        return EXIT_FAILURE;
    }
    if (argc > 4 && parse_int(argv[4], "nsteps", &nsteps) != 0) {
        return EXIT_FAILURE;
    }
    if (nsteps < 2) {
        fprintf(stderr, "nsteps must be at least 2\n");
        return EXIT_FAILURE;
    }

    if (tov_load_eos(eos_file, &eos) != 0) {
        return EXIT_FAILURE;
    }

    rho_h = (rho_end - rho_start) / (double)(nsteps - 1);
    rho_tmp = rho_start;

    tov_print_header(stdout);

    /* Loop over central density */
    for (int i = 0; i < nsteps; i++) {
        TovResult result;
        double rhoc = rho_tmp * TOV_FM3CM3;

        if (TOV_SOLVER(&eos, rhoc, &result) != 0) {
            tov_free_eos(&eos);
            return EXIT_FAILURE;
        }
        tov_print_result(stdout, &result);
        rho_tmp = rho_tmp + rho_h;
    }

    tov_free_eos(&eos);
    return EXIT_SUCCESS;
}
