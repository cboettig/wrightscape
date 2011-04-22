/**
 * @file optimizers.h
 * @author Carl Boettiger <cboettig@gmail.com>
 * @date 22 April 2011
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_siman.h>

/* set up parameters for this simulated annealing run */
#define N_TRIES 200     	/* how many points do we try before stepping */
#define ITERS_FIXED_T 200	/* how many iterations for each T? */
#define STEP_SIZE .02 		/* max step size in random walk */
#define K 1.0				/* Boltzmann constant */
#define T_INITIAL 0.05		/* initial temperature */
#define MU_T 1.004		    /* damping factor for temperature */
#define T_MIN .008

/* setup parameters for multimin method*/
#define INIT_STEP .2
#define MAX_ITER 10000
#define ERR_TOL 1e-8
#define PRINT 1
#include <gsl/gsl_multimin.h>


double optim_func (const gsl_vector *v, void *params);
double multimin(gsl_vector *x, void * params);
double siman(gsl_vector * x, void * params, gsl_rng * rng); 

