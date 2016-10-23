//
//  Copyright © 2016 Lennart Oymanns. All rights reserved.
//

#ifndef vegas_h
#define vegas_h

struct vegas_state;

/**
 vegas_new creates a vegas_state object. The parameter ndim is the 
 dimension of the integration and nbin is the number of bins per 
 dimension (try e.g. 50). The memory of the returned 
 vegas_state object has to be freed with the function vegas_free().
 */
struct vegas_state *  vegas_new(int ndim, int nbin);

/**
 vegas_free frees the memory of vegas_state object.
 */
void vegas_free(struct vegas_state * state);

/**
 vegas_add adds a integrand value to the vegas state.
 */
void vegas_add(struct vegas_state * state,
               double integrand,
               double wgt,
               int use);

/**
 vegas_get_integrand_args returns integrand arguments. 
 The array u has to contain ndim random numbers uniformly distributed in [0,1).
 The random numbers u are transformed and are then written to the output array x.
 The vegas weight is written to weight. The values in x are then to be used as 
 input variables for the integrant.
 vegas_get_integrand_args returns 1, if there are still iterations to perform.
 It returns 0 if the integration is finished. weight and x are not changed in 
 that case.
 */
int vegas_get_integrand_args(struct vegas_state * state,
                             const double u[],
                             double *weight,
                             double * x);

/**
 vegas_get_integral sets the variables integral, error and chi2_ndf to the
 integral, the error and to chi^2/ndf.
 */
void vegas_get_integral(struct vegas_state * state,
                        double * integral,
                        double * error,
                        double * chi2_ndf);

/**
 vegas_init_integration resets state to a clean state without any integrand
 evaluation. The grid is unchanged. The number of iterations and points per 
 iteration are set by the parameters interations and nperit.
 */
void vegas_init_integration(struct vegas_state * state,
                            int iterations,
                            int nperit);

void vegas_set_verbose(struct vegas_state * state, int lvl);

#endif /* vegas_h */
