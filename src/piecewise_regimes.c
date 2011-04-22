/** 
 * @file piecewise_regimes.c
 * @author Carl Boettiger <cboettig@gmail.com>
 * @date 17 May 2010
 *
 * @section LICENSE 
 *
 * GPLv3
 *
 * @section DESCRIPTION
 *
 * Code for computing the likelihood of a multitype OU process on a phylogenetic tree
 *
 */

#include "likelihood.h"

/**
* @brief Function which wraps general format used by siman and multimin 
*        into the likelihood format of calc_lik
*
* @param v vector of values being optimized by the algorithm 
* @param params specifies the tree and observed data
*
* @return minus log likelihood
*/
double optim_func (const gsl_vector *v, void *params)
{
	tree * mytree = (tree *) params;
	int i, n_regimes = mytree->n_regimes;
// Force root to match Theta of regime assigned to it, should also reduce vector size!
	mytree->Xo[0] = mytree->theta[mytree->regimes[0]];
 //	mytree->Xo[0] = gsl_vector_get(v, 0);
	for(i = 0; i < n_regimes; i++){
		mytree->alpha[i] = v->data[i];
		mytree->theta[i] = v->data[n_regimes+i];
		mytree->sigma[i] = v->data[2*n_regimes+i];

		if( (mytree->alpha[i] < 0) | (mytree->sigma[i] < 0) ){ return GSL_POSINF; }
		if(  mytree->theta[i] < 0 ){ return GSL_POSINF; }
	}
	
	return -calc_lik(mytree->Xo, mytree->alpha, mytree->theta, mytree->sigma,
    			     mytree->regimes, mytree->ancestor, mytree->branch_length,
                     mytree->traits, &(mytree->n_nodes), mytree->lca_matrix); 
}


/** Fit the model by maximum likelihood, called by R */
void fit_model(double * Xo, 
	double * alpha, 
	double * theta, 
	double * sigma, 
	int * regimes, 
	int * ancestor, 
	double * branch_length, 
	double * traits, 
	int * n_nodes, 
	int * n_regimes,
	double * llik,
	int * use_siman)
{

//  Perhaps don't want to see gsl_errors once we're handling in R?
//	gsl_set_error_handler_off ();
	
    /** First, initialize the tree structure to facilitate passing all these elements */
	int i,j;
	tree * mytree = (tree  *) malloc(sizeof(tree));
	mytree->Xo = Xo;
	mytree->alpha = alpha;
	mytree->theta = theta;
	mytree->sigma = sigma;
	mytree->regimes = regimes;
	mytree->ancestor = ancestor;
	mytree->branch_length = branch_length;
	mytree->traits = traits;
	mytree->n_nodes = *n_nodes;
	mytree->n_regimes = *n_regimes;

    /** Second, save time by storing the least common ancestor matrix */
	/* Save time by calculating least common ancestor ahead of time.
	 * Though this does the calc for all nodes, only tip ones are used. 
	 * Current get_lca algorithm is only designed for tips anyway.  */
	double sep; // separation time, not actually used
	mytree->lca_matrix = (int *) malloc(gsl_pow_2(*n_nodes) * sizeof(int) );
	for(i=0; i < *n_nodes; i++){
		for(j=0; j < *n_nodes; j++){
			mytree->lca_matrix[*n_nodes*i+j] = get_lca(i,j, *n_nodes, ancestor, branch_length, &sep);
		}
	}

    /** Then initialize the vector of parameters to be searched */
	gsl_vector *x = gsl_vector_alloc(3 * *n_regimes);
//	gsl_vector_set(x, 0, *Xo);  /* Forces Xo to match theta (?) */
	for(i = 0; i < *n_regimes; i++){
		gsl_vector_set(x, i, alpha[i]);
		gsl_vector_set(x, *n_regimes+i, theta[i]);
		gsl_vector_set(x, 2 * *n_regimes+i, mytree->sigma[i]);
	}

    /** Select and call a likelihood optimization routine 
     *  given the tree (with data) and parameters */
	gsl_rng * rng = gsl_rng_alloc(gsl_rng_default);
	if(*use_siman)
		*llik = siman(x, mytree, rng); 
	else 
		*llik = multimin(x, mytree);


	gsl_rng_free(rng);
	gsl_vector_free(x);
	free(mytree->lca_matrix);
	free(mytree);
}




/* Externally exported function to determine LCA matrix */
/**
* @brief Determine the last common ancestor matrix of the tips, used in the likelihood calculation
*
* @param ancestor ancestor of the specified nide
* @param branch_length double branch lengths below each of the nodes
* @param n_nodes: integer number of nodes (including tips) in the tree 
* @param lca_matrix must be an integer array of size n_nodes^2
*/
void calc_lca(int * ancestor, double * branch_length, 
              int * n_nodes, int * lca_matrix)
{
	int i,j;
	/* Save time by calculating least common ancestor ahead of time.
	 * Though this does the calc for all nodes, only tip ones are used. 
	 * Current get_lca algorithm is only designed for tips anyway.  */
	double sep; // separation time, not actually used but lca_matrix wants a handle to it.  
	lca_matrix = (int *) malloc(gsl_pow_2(*n_nodes) * sizeof(int) );
	for(i=0; i < *n_nodes; i++){
		for(j=0; j < *n_nodes; j++){
			lca_matrix[*n_nodes*i+j] = get_lca(i,j, *n_nodes, ancestor, branch_length, &sep);
		}
	}

}




/** Simulate under the model.  Called by R */
void simulate_model (double * Xo, double * alpha, double * theta, 
                     double * sigma, int * regimes, int * ancestor, 
                     double * branch_length, double * traits, int * n_nodes,
                     int * n_regimes, double * llik, double * seed)
{

	int i,j;
	tree * mytree = (tree  *) malloc(sizeof(tree));
	mytree->Xo = Xo;
	mytree->alpha = alpha;
	mytree->theta = theta;
	mytree->sigma = sigma;
	mytree->regimes = regimes;
	mytree->ancestor = ancestor;
	mytree->branch_length = branch_length;
	mytree->traits = traits;
	mytree->n_nodes = *n_nodes;
	mytree->n_regimes = *n_regimes;

	/* Save time by calculating least common ancestor ahead of time. */
	double sep; // separation time, not actually used
	mytree->lca_matrix = (int *) malloc(gsl_pow_2(*n_nodes) * sizeof(int) );
	for(i=0; i < *n_nodes; i++){
		for(j=0; j < *n_nodes; j++){
			mytree->lca_matrix[*n_nodes*i+j] = get_lca(i,j, *n_nodes, ancestor, branch_length, &sep);
		}
	}


    /* Allocate the parameter vector */
	gsl_vector *x = gsl_vector_alloc(1 + 3 * *n_regimes);
	gsl_vector_set(x, 0, *Xo);
	for(i = 0; i < *n_regimes; i++){
		gsl_vector_set(x, 1+i, alpha[i]);
		gsl_vector_set(x, 1+*n_regimes+i, theta[i]);
		gsl_vector_set(x, 1+2 * *n_regimes+i, mytree->sigma[i]);
	}
	
	gsl_rng * rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, *seed);

   /* These are the only lines that differ from the fit_models wrapper */ 
	simulate (rng, mytree);
	*llik =  calc_lik(mytree->Xo, mytree->alpha, mytree->theta, mytree->sigma, 
			          mytree->regimes, mytree->ancestor, mytree->branch_length, 
                      mytree->traits, &(mytree->n_nodes), mytree->lca_matrix);

//	for(i=0; i< mytree->n_regimes; i++) printf("%g %g %g\n", mytree->alpha[i], mytree->theta[i], mytree->sigma[i] );
//	printf("sim llik: %g\n", *llik);


	gsl_rng_free(rng);
	gsl_vector_free(x);
	free(mytree->lca_matrix);
	free(mytree);
}





/** Test code by running in pure C */
int main(void)
{
/*	
	double branch_length[]	= {0, 1, 1, 3, 2, 1, 1}; 
	const int ancestor[]	={-1, 2, 0, 0, 2, 1, 1};
	const double traits[]	= {0, 0, 0, 5, 2, 1, 3};	
	const int regimes[]		= {0, 0, 0, 1, 0, 0, 0};
	int n_nodes = 7;

	double Xo = 1;
	double alpha[2] = {2, 2};
	double theta[2] = {1, 5}; 
	double sigma[2] = {1, 1};
/ */
	int n_nodes = 45;
	double branch_length[] = { 0.0, 12./38, 20./38,  2./38,  2./38,  4./38,  
							 8./38,  5./38,  5./38,  5./38,  5./38, 10./38,  
							 9./38,  4./38,  8./38,  2./38, 20./38,  2./38,  
							 4./38,  2./38, 1./38,  2./38, 26./38, 4./38, 
							 2./38, 2./38,  2./38,  2./38, 15./38, 10./38, 
							 10./38, 10./38, 10./38, 16./38, 12./38,  4./38,  
							 2./38,  2./38, 10./38,  8./38,  2./38,  1./38,  
							 1./38,  2./38,  2./38};

	int ancestor[] = {-1,  0,  1,  2,  3,  2,  0,  6,  7,  8,  9,
						8,  7, 12, 13,14,  6, 16, 17, 18, 19, 18, 1,
						3,  4,  4,  5,  5,  9, 10, 10, 11, 11, 12, 
						13, 14, 15, 15, 16, 17, 19, 20, 20, 21, 21};

	double traits[] = {   0.0000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
						  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
						  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.602690, 2.660260,
						  2.660260, 2.653242, 2.674149, 2.701361, 3.161247, 3.299534, 3.328627, 3.353407,
						  3.360375, 3.049273, 2.906901, 2.980619, 2.933857, 2.975530, 3.104587, 3.346389,
						  2.928524, 2.939162, 2.990720, 3.058707, 3.068053};
	//int regimes[45] = {0};
	int regimes[] =  {1, 1, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
						2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1}; 

/*
	// ou1 parameters
	double Xo = 2.953806;
	double alpha[3] = { 0.192155, 0.192155, 0.192155};
	double theta[3] = {2.953806,  2.953806, 2.953806}; 
	double sigma[3] = {sqrt( .048365), sqrt( .048365), sqrt( .048365)};

	double Xo = 1;
	double alpha[3] = {1., 1., 1.};
	double theta[3] = {3., 3., 3.};
	double sigma[3] = {1.,  1., 1. };

*/

	double Xo = 3;
	double alpha[3] = {2.6, 2.6, 2.6};
	double theta[3] = {3., 3.,3};
	double sigma[3] = {.2,  .2, .2 };
	int n_regimes = 3;
	double llik = 0;

	int use_siman = 0;

	fit_model(&Xo, alpha, theta, sigma, regimes, ancestor, branch_length, traits, &n_nodes, &n_regimes, &llik, &use_siman);
	printf("Xo = %g\n", Xo);
	printf("alphas: %g %g %g\n", alpha[0], alpha[1], alpha[2]);
	printf("thetas: %g %g %g\n", theta[0], theta[1], theta[2]);
	printf("sigmas: %g %g %g\n", sigma[0], sigma[1], sigma[2]);

	printf("log likelihood: %g\n", llik);

	double seed = 1.0;
	printf("seed: %lo\n", (unsigned long int) seed);

	simulate_model(&Xo, alpha, theta, sigma, regimes, ancestor, branch_length, traits, &n_nodes, &n_regimes, &llik, &seed);

//	int i;
//	for(i=0; i < 45; i++) printf("%g\n", traits[i]);

	
	return 0;
}
