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
void traitmodels(	const double * times, ///< The lengths of each node. Nodes identified by index value.  
					const int * ancestors, ///< A list of the ancestors of each node.  The root node is number 0.  
					const double * traits, ///< The trait data
					const int * states, ///< Regimes list for each node
					const int * nstates, ///< number of regimes (currently must be 1)
					const int * n, ///< number of nodes
					const double * pars, ///< list of parameters.  Sigma, alpha, theta.  
					const int * fitpars, ///< logical indicating which parameters are fit and which held constant
					const int * npars, ///< how many parameters total
					double * likelihood
				)
{
	int i;
	tree * mytree = tree_alloc(*n, *npars);
	tree_init(mytree, times, ancestors, traits, states, *nstates, *n, pars, fitpars, GRIDSIZE);

//	tree_print(mytree);
//	matrix_regimes(mytree);
	
//	printf("analytic bm llik = %lf\n", bm_likelihood(mytree) );
//	printf("analytic ou llik = %lf\n", ou_likelihood(mytree) );
//	printf("matrix llik = %lf\n",  matrix_likelihood(mytree) );


	/* consider using a pointer to the likelihood method 
	 * siman method needs more testing */
/*	gsl_rng * rng = gsl_rng_alloc(gsl_rng_default);
	siman(mytree, rng);
	gsl_rng_free(rng);
*/
	*likelihood = multimin(mytree);
	tree_free(mytree);
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
/*	// test data
	const double times[] =  {0, 1, 1, 3, 2, 1, 1}; 
	const int ancestors[] = {-1, 2, 0, 0, 2, 1, 1};
	const double traits[] = {0, 0, 0, -2, 1, 2, 2};	
	const int states[] =	{0, 0, 0, 0, 0, 0, 0};
	const int nstates = 1;
	const int n = 7;
*/

	// sigma2 alpha, theta
	const double pars[] = {.04836469, 0.1921554, 2.953806}; //, 4.1, .21, 3.1, 0.1, 0.1}; 
	// Which parameters should be fitted??
	const int fitpars[] = {1,   1,   0}; //,  0,   0,   0,   0,   0};
	const int npars = 3;
	double likelihood = 0;

//	traitmodels(times, ancestors, traits, states, &nstates, &n, pars, fitpars, &npars, &likelihood);

	// bimac data
	const int bimac_n = 45;
//	const double bimac_times[] = { 0, 12, 20,  2,  2,  4,  8,  5,  5,  5,  5, 10,  9,  4,  8,  2, 20,  2,  4,  2, 1,  2, 26,
//							 4, 2, 2,  2,  2, 15, 10, 10, 10, 10, 16, 12,  4,  2,  2, 10,  8,  2,  1,  1,  2,  2};
	double bimac_times[] = { 0, 12./38, 20./38,  2./38,  2./38,  4./38,  8./38,  5./38,  5./38,  5./38,  5./38, 10./38,  9./38,  4./38,  8./38,  2./38, 20/38,  2./38,  4./38,  2./38, 1./38,  2./38, 26./38, 4./38, 2./38, 2./38,  2./38,  2./38, 15./38, 10./38, 10./38, 10./38, 10./38, 16./38, 12./38,  4./38,  2./38,  2./38, 10./38,  8./38,  2./38,  1./38,  1./38,  2./38,  2./38};

//	gsl_vector_view A = gsl_vector_view_array(bimac_times, 45);
//	gsl_vector_fprintf(stdout, &A.vector, "%g");

	const int bimac_ancestors[] = {-1,  0,  1,  2,  3,  2,  0,  6,  7,  8,  9,  8,  7, 12, 13,14,  6, 16, 17, 18, 19, 18, 1,
							3,  4,  4,  5,  5,  9, 10, 10, 11, 11, 12, 13, 14, 15, 15, 16, 17, 19, 20, 20, 21, 21};

	const double bimac_traits[] = { 2.900000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
							  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
							  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.602690, 2.660260,
							  2.660260, 2.653242, 2.674149, 2.701361, 3.161247, 3.299534, 3.328627, 3.353407,
							  3.360375, 3.049273, 2.906901, 2.980619, 2.933857, 2.975530, 3.104587, 3.346389,
							  2.928524, 2.939162, 2.990720, 3.058707, 3.068053};
	const int bimac_states[45] = {0};
	const int bimac_nstates = 1;

	const int LP_states[] =  {1, 1, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
							2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1}; 

	traitmodels(bimac_times, bimac_ancestors, bimac_traits, bimac_states, &bimac_nstates, &bimac_n, pars, fitpars, &npars, &likelihood);

	return 0;
}


