/**
 * @file jointprob.c 
 * @author Carl Boettiger, <cboettig@gmail.com>
 *
 * @section DESCRIPTION
 *  
 *
 * @section LICENSE
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 */

#include "tree.h"
/** 
 * The function traitmodels() will be called directly from R 
 * (1) assembles the data from R into the appropriate data structure by calling the tree function
 * (2) initializes necessary things like random number generator
 * (3) calls the parameter search algorithm (minimizer or MCMC)
 * (4) cleans up and returns data
 */
void traitmodels(const double * times, const int * ancestors, const double * traits, const int * states, const int * nstates, const int * n, const double * pars, const int * fitpars, const int * npars)
{
	int i;
	tree * mytree = tree_alloc(*n, *npars);
	tree_init(mytree, times, ancestors, traits, states, *nstates, *n, pars, fitpars, GRIDSIZE);

//	tree_print(mytree);
//	matrix_regimes(mytree);

	int nfreepars = 0; // DoF for AIC scores
	for(i=0;i<*npars;i++) nfreepars += fitpars[i];
	
	printf("analytic bm llik = %lf\n", bm_likelihood(mytree) );
	printf("analytic ou llik = %lf\n", ou_likelihood(mytree) );
	printf("matrix llik = %lf\n",  matrix_likelihood(mytree) );

	gsl_rng * rng = gsl_rng_alloc(gsl_rng_default);

	/* consider using a pointer to the likelihood method */
//	siman(mytree, rng);
//	multimin(mytree);

	tree_free(mytree);
	gsl_rng_free(rng);
}

	

/** 
 * The main function simply provides the data to traitmodels().  
 * Eventually this data will instead be handed from the R wrapper to trait models.  
 * 
 * Two datasets are currently provided here.  A toy example of 4 tips for testing the methods and
 * the bimac Anoles data provided in ouch R program (Butler & King 2004).  
 *
 * The data is specified by a tree (ancestors list determines topology) and branch lengths (times),
 * together with a list of states at each node in the tree.  Traits at internal nodes aren't actually used.  
 * The root trait will be treated as a parameter.  
 *
 *
 * */
int main(void)
{
	// test data
	const double times[] =  {0, 1, 1, 3, 2, 1, 1}; 
	const int ancestors[] = {-1, 2, 0, 0, 2, 1, 1};
	const double traits[] = {0, 0, 0, -2, 1, 2, 2};	
	const int states[] =	{0, 0, 0, 0, 0, 0, 0};
	const int nstates = 1;
	const int n = 7;


	const double pars[] = {4., .02, 2.5}; //, 4.1, .21, 3.1, 0.1, 0.1}; 
	const int fitpars[] = {1,   1,   1}; //,  0,   0,   0,   0,   0};
	const int npars = 3;


	traitmodels(times, ancestors, traits, states, &nstates, &n, pars, fitpars, &npars);

	// bimac data
	const int bimac_n = 45;
	const double bimac_times[] = { 0, 12, 20,  2,  2,  4,  8,  5,  5,  5,  5, 10,  9,  4,  8,  2, 20,  2,  4,  2, 1,  2, 26,
							 4, 2, 2,  2,  2, 15, 10, 10, 10, 10, 16, 12,  4,  2,  2, 10,  8,  2,  1,  1,  2,  2};
	const int bimac_ancestors[] = {-1,  0,  1,  2,  3,  2,  0,  6,  7,  8,  9,  8,  7, 12, 13,14,  6, 16, 17, 18, 19, 18, 1,
							3,  4,  4,  5,  5,  9, 10, 10, 11, 11, 12, 13, 14, 15, 15, 16, 17, 19, 20, 20, 21, 21};

	const double bimac_traits[] = { 2.400000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
							  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
							  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.602690, 2.660260,
							  2.660260, 2.653242, 2.674149, 2.701361, 3.161247, 3.299534, 3.328627, 3.353407,
							  3.360375, 3.049273, 2.906901, 2.980619, 2.933857, 2.975530, 3.104587, 3.346389,
							  2.928524, 2.939162, 2.990720, 3.058707, 3.068053};
	const int bimac_states[45] = {0};
	const int bimac_nstates = 1;

	const int LP_states[] =  {1, 1, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
							2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1}; 

	traitmodels(bimac_times, bimac_ancestors, bimac_traits, bimac_states, &bimac_nstates, &bimac_n, pars, fitpars, &npars);

	return 0;
}


