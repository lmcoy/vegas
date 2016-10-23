//
//  main.c
//  vegas
//
//  Created by Lennart Oymanns on 19/06/16.
//  Copyright © 2016 Lennart Oymanns. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <mpi.h>

#include "vegas.h"

void get_random(int n, double * u) {
    int i = 0;
    for (i = 0; i < n; i++) {
        u[i] = drand48();
    }
}

void f() {
    struct vegas_state * state = vegas_new(2, 50);
    vegas_set_verbose(state, 1);
    
    double x[2];
    double r[2];
    get_random(2, r);
    vegas_init_integration(state, 5, 10000);
    
    double wgt = 1.0;
    while(vegas_get_integrand_args(state, r, &wgt, x)) {
        double f = x[0] * x[1];
        vegas_add(state, f, wgt, 1);
        get_random(2, r);
    }
    
    vegas_init_integration(state, 5, 100000);
    while(vegas_get_integrand_args(state, r, &wgt, x)) {
        double f = x[0] * x[1];
        vegas_add(state, f, wgt, 1);
        get_random(2, r);
    }
    double integral = 0.0;
    double error = 0.0;
    double chi2 = 0.0;
    vegas_get_integral(state, &integral, &error, &chi2);
    
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf ("deviation in sigma: %g\n", (integral - 0.25)/error);
    }
    vegas_free(state);
}

int main(int argc, const char * argv[]) {
    MPI_Init(NULL,NULL);
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand48(time(NULL) + 13 * rank);
    f();
    
    MPI_Finalize();

    return 0;
}
