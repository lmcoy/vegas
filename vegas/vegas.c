//
//  Copyright © 2016 Lennart Oymanns. All rights reserved.
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

//#define NOMPI

#ifndef NOMPI
#include <mpi.h>
#endif

#include "vegas.h"
#include "matrix.h"

static int min(int a, int b) {
    if (a < b) {
        return a;
    }
    return b;
}

struct int_state {
    int *ia;
    struct Matrix * d;
    double fb;
    double f2b;
    int ncalls;
};

static void int_state_reset(struct int_state *istate) {
    matrix_set_zero(istate->d);
    istate->fb = 0.0;
    istate->f2b = 0.0;
    istate->ncalls = 0;
}

static void int_state_init(struct int_state *istate, int ndim, int nbin) {
    istate->d = matrix_new(nbin, ndim);
    istate->ia = malloc(sizeof(int) * ndim);
    
    int_state_reset(istate);
}

static void int_state_free(struct int_state *istate) {
    free(istate->ia);
    istate->ia = NULL;
    matrix_free(istate->d);
}

static void rebin(const double rc, const int nd, const double *r,
                  struct Matrix *xi, const int j);

void vegas_reset_grid(struct vegas_state *state);
void vegas_init_grid(struct vegas_state *state);

struct vegas_state {
    struct int_state Istate;
    int Iterations; /// number of iterations
    int NperIt; /// points per iteration
    int Ndim; /// integral dimension
    int Nbin; /// number of bins per dimension
    struct Matrix *xi;
    struct {
        int n; /// current point
        int it; /// current iteration
    } State;
    struct {
        int rank; /// MPI rank
        int nproc; /// number of processes
    } MPI;
    double si; /// sum of weighted iteration results
    double schi;
    double swgt; /// sum of iteraton weights
    int refine_grid;
    int verbose;
    clock_t start_time;
};

void refine_grid(struct vegas_state * state, struct Matrix *d);

struct vegas_state *  vegas_new(int ndim, int nbin) {
    struct vegas_state * state = malloc(sizeof(struct vegas_state));
    state->Iterations = 0;
    state->NperIt = 0;
    
    
    state->xi = matrix_new(ndim, nbin);
#ifndef NOMPI
    MPI_Comm_rank(MPI_COMM_WORLD, &state->MPI.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &state->MPI.nproc);
#else
    state->MPI.rank = 0;
    state->MPI.nproc = 1;
#endif
    state->Ndim = ndim;
    state->Nbin = nbin;
    
    state->si = 0.0;
    state->schi = 0.0;
    state->swgt = 0.0;
    
    state->refine_grid = 1;
    state->State.n = 0;
    state->State.it = 0;
    state->verbose = 0;
    state->start_time = 0;
    
    int_state_init(&state->Istate, ndim, nbin);
    
    vegas_reset_grid(state);
    vegas_init_grid(state);
    
    return state;
}

void vegas_reset_grid(struct vegas_state *state) {
    int ndim = state->Ndim;
    int j;
    for (j = 0; j < ndim; j++) {
        matrix_set(state->xi, j, 0, 1.0);
    }
}

void vegas_init_grid(struct vegas_state *state) {
    int ndim = state->Ndim;
    int nbin = state->Nbin;
    double *r = malloc(sizeof(double) * nbin);
    int i, j;
    for (i = 0; i < nbin; i++) {
        r[i] = 1.0;
    }
    for (j = 0; j < ndim; j++) {
        rebin(1.0 / (double)nbin, nbin, r, state->xi, j);
    }
    free(r);
    r = NULL;
}

void vegas_free(struct vegas_state * state) {
    matrix_free(state->xi);
    state->xi = NULL;
    int_state_free(&state->Istate);
    free(state);
    state = NULL;
}

int vegas_get_integrand_args(struct vegas_state * state, const double u[], double *weight, double * x) {
    if (state->State.it>= state->Iterations) {
        return 0;
    }
    double wgt = 1.0;

    int ndim = state->Ndim;
    int nbin = state->Nbin;
    struct int_state * istate = &state->Istate;
    int j = 0;
    for (j = 0; j < ndim; j++) {
        double xn = (1 - u[j]) * nbin + 1.0;
        istate->ia[j] = min((int)(xn), nbin);
        double x_1 = matrix_get(state->xi, j, istate->ia[j] - 1);
        double x_2 = 0.0;
        if (istate->ia[j] > 1) {
            x_2 = matrix_get(state->xi, j, istate->ia[j] - 2);
        }
        double delta_x = x_1 - x_2;
        x[j] = x_2 + (xn - istate->ia[j]) * delta_x;
        wgt *= delta_x * (double)nbin;
    }
    *weight = wgt;
    return 1;
}

static void finish_iteration(struct vegas_state * state, int update_grid) {
    static const double TINY = 1.0e-30;
    struct int_state * istate = &state->Istate;

    state->State.n = 0;
    state->State.it += 1;
    int nbin = state->Nbin;
    int ndim = state->Ndim;
    
#ifndef NOMPI
    double fb = 0.0;
    MPI_Reduce(&istate->fb, &fb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double f2b = 0.0;
    MPI_Reduce(&istate->f2b, &f2b, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    struct Matrix *d = NULL;
    if (update_grid) {
        d = matrix_new(nbin, ndim);
        MPI_Reduce(istate->d->m, d->m, nbin*ndim, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    }
    
    int ncalls_total = 0;
    MPI_Reduce(&istate->ncalls, &ncalls_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    double fb = istate->fb;
    double f2b = istate->f2b;
    struct Matrix * d = istate->d;
    int ncalls_total = istate->ncalls;
#endif
    
    if(state->MPI.rank == 0) {
        double jac = 1.0/((double)ncalls_total);
        double dv2g = 1.0 / (ncalls_total - 1.0);
        fb *= jac;
        f2b *= jac * jac;
        if (update_grid) {
            matrix_mul_factor(d, jac * jac);
        }
        f2b = sqrt(f2b * ncalls_total);
        f2b = (f2b - fb) * (f2b + fb);
        if (f2b <= 0.0) {
            f2b = TINY;
        }
        double integral_it = fb;
        double sigma_it = f2b * dv2g;
        double wgt = 1.0 / sigma_it;
        state->si += wgt * integral_it;
        state->schi += wgt * integral_it * integral_it;
        state->swgt += wgt;
        double integral = state->si / state->swgt;
        int it = state->State.it;
        double chi2_ndf = (state->schi - state->si * integral) / (it + 0.0001);
        if (chi2_ndf < 0.0) {
            chi2_ndf = 0.0;
        }
        double sd = sqrt(1.0 / state->swgt);
        sigma_it = sqrt(sigma_it);
        
        if (state->verbose) {
            clock_t end = clock();
            double elapsed_secs = (double)(end - state->start_time) / CLOCKS_PER_SEC;
            time_t now = time(0);
            int it_left = state->Iterations - it;
            now += it_left * elapsed_secs;
            char time_buffer[30];
            struct tm now_tm = *localtime(&now);
            strftime(time_buffer, 30, "%c", &now_tm);
            
            state->start_time = end;
            printf( "iteration %d:\n", it);
            printf( "  %g +- %g (chi^2/ndf = %g) %s\n", integral, sd, chi2_ndf, time_buffer);
        }
        if (update_grid) {
            refine_grid(state, d);
        }
    }
#ifndef NOMPI
    // broadcast new grid
    if (update_grid) {
        matrix_free(d);
        MPI_Bcast(state->xi->m, ndim * nbin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&state->si, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&state->schi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&state->swgt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void vegas_add(struct vegas_state * state, double integrand, double wgt, int use) {
    if (state->start_time == 0) {
        state->start_time = clock();
    }
    double f = wgt * integrand;
    double f2 = f * f;
    struct int_state * istate = &state->Istate;
    istate->fb += f;
    istate->f2b += f2;
    istate->ncalls += 1;
    int ndim = state->Ndim;
    int j = 0;
    for (j = 0; j < ndim; j++) {
        int iaj = istate->ia[j];
        matrix_add_to_elem(istate->d, iaj - 1, j, f2);
    }
    if (!use) {
        return;
    }
    
    state->State.n += 1;
    if( state->State.n >= state->NperIt ) {
        finish_iteration(state, 1);
    }

}

void vegas_init_integration(struct vegas_state * state, int iterations, int nperit) {
    state->Iterations = iterations;
    
    int nproc = state->MPI.nproc;
    int rank = state->MPI.rank;
    state->NperIt = nperit/nproc;
    int rem = nperit % nproc;
    if (rem>0 && rank <= rem) {
        state->NperIt += 1;
    }
    state->State.it = 0;
    state->State.n = 0;
    state->si = 0.0;
    state->schi = 0.0;
    state->swgt = 0.0;
    state->start_time = 0;
}

void vegas_get_integral(struct vegas_state * state, double * integral, double * error, double *chi2_ndf) {
    *integral = state->si / state->swgt;
    int it = state->State.it;
    *chi2_ndf = (state->schi - state->si * *integral) / (it + 0.0001);
    if (*chi2_ndf < 0.0) {
        *chi2_ndf = 0.0;
    }
    *error = sqrt(1.0 / state->swgt);
}

void vegas_set_verbose(struct vegas_state * state, int lvl) {
    state->verbose = lvl;
}

static void rebin(const double rc, const int nd, const double *r,
                  struct Matrix *xi, const int j) {
    int k = 0;
    double dr = 0.0, xn = 0.0, xo = 0.0;
    
    double * xin = malloc(sizeof(double)*(nd-1));
    
    int i;
    for (i = 0; i < nd - 1; i++) {
        while (rc > dr) {
            dr += r[k];
            k += 1;
        }
        if (k > 1) {
            xo = matrix_get(xi, j, k - 2);
        }
        xn = matrix_get(xi, j, k - 1);
        dr -= rc;
        xin[i] = xn - (xn - xo) * dr / r[k - 1];
    }
    for (i = 0; i < nd - 1; i++) {
        matrix_set(xi, j, i, xin[i]);
    }
    free(xin);
    xin = NULL;
    matrix_set(xi, j, nd - 1, 1.0);
}

void refine_grid(struct vegas_state * state, struct Matrix *d) {
    static const double ALPH = 1.5, TINY = 1.0e-30;
    if (!state->refine_grid) { return; }
/* refine grid */
    int ndim = state->Ndim;
    int nbin = state->Nbin;
    double * dt = malloc(sizeof(double)*ndim);
    double * r = malloc(sizeof(double)*nbin);
    int i = 0;
    int j;
    for (j = 0; j < ndim; j++) {
        double xo = matrix_get(d, 0, j);
        double xn = matrix_get(d, 1, j);
        matrix_set(d, 0, j, (xo + xn) / 2.0);
        dt[j] = matrix_get(d, 0, j);
        for (i = 2; i < nbin; i++) {
            double rc = xo + xn;
            xo = xn;
            xn = matrix_get(d, i, j);
            matrix_set(d, i - 1, j, (rc + xn) / 3.0);
            dt[j] += matrix_get(d, i - 1, j);
        }
        matrix_set(d, nbin - 1, j, (xo + xn) / 2.0);
        dt[j] += matrix_get(d, nbin - 1, j);
    }
    for (j = 0; j < ndim; j++) {
        double rc = 0.0;
        for (i = 0; i < nbin; i++) {
            double dij = matrix_get(d, i, j);
            if (dij < TINY) {
                matrix_set(d, i, j, TINY);
                dij = TINY;
            }
            r[i] = pow((1.0 - dij / dt[j]) / (log(dt[j]) - log(dij)),
                       ALPH);
            rc += r[i];
        }
        rebin(rc / (double)nbin, nbin, r, state->xi, j);
    }
    free(dt);
    dt = NULL;
    free(r);
    r = NULL;
}


