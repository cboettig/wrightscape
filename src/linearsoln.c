/** 
 * @file linearsoln.c
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
 *  
 *  Note: Comments formatted for Doxygen, which generates the documentation files.  
 *
  */

#include "tree.h"


/** The branch length between nodes i and j; (twice their divergence time; down to common ancestor and back).  
 *  Currently this requires i & j to be the same age (i.e. both are tips in an ultrametric tree).  The code 
 *  opens with a warning message to handle this.  
 *
 *  This should really implement a clever and general algorithm.  See the wikipedia page on Least Common Ancestor
 *  or the tutorial, \see{http://www.topcoder.com/tc?module=Static&d1=tutorials&d2=lowestCommonAncestor}
 *
 *  Really no need to be fast as it should only be called once per tree, as it is indep of parameters and painting!
 * */
double sep_time(int i, int j, tree * mytree)
	{
	if(!mytree->tip[i] || !mytree->tip[j]) { printf("currently i and j must be tips!\n"); return 0; }

	int n = mytree->n;
	double * time_to_ancestor = (double *) malloc(n * sizeof(double));
	int * ancestor_list = (int *) malloc(n * sizeof(int));
	int k = 0, s = 0;
	time_to_ancestor[0] = 0.0;
	double sep;

	for(k=0;k<n;k++) ancestor_list[k]=0;
	
	k = 0;
	while(1){
		ancestor_list[k] = i;
		time_to_ancestor[k+1] = time_to_ancestor[k] + mytree->time[i];
		if(i==0) break;
		i = mytree->ancestor[i];
		k++;
	}
	while(1){
		for(k=0;k<n;k++){
			if(j == ancestor_list[k]){ 
				s = j;
				j = 0;
				break; 
			}
		}
		if(j==0) break;
		j = mytree->ancestor[j];
	}
//	printf("LCA = %d, time = %lf \n", s, time_to_ancestor[k] );
	
	sep = time_to_ancestor[k];	
	free(time_to_ancestor);
	free(ancestor_list);
	return 2*sep;
}



/**
 * 
 * The likelihood of a model under the OU process,
 * @f[ dX_t = \alpha (\theta - X_t) \d t + \sigma \d W_t @f]
 * is multivariate normal. If only the n tips are specified (not computing likelihood of ancestral state) then
 * this is determined by the n mean values and the n x n covariance matrix.  The expected value at any node is 
 * determined only by the age of the node and the intitial condition of the root node (@f$ X_0 @f$):
 * @f[ E(X_t | X_0 ) = X_0 e^{-\alpha t} + \theta (1- e^{-\alpha t}) @f] 
 * The i,j element of the covariance matrix is
 * @f[ V_{ij} = \frac{\sigma^2}{2\alpha} (1 - e^{-2 \alpha s_{ij} }) e^{-2\alpha (t-s_{ij} )} @f]
 * where the nodes are all measured at time t (present day) and have diverged for a time \f$ t-s_{ij} \f$ 
 * (thus shared time \f$ s_{ij} \f$ since \f$ X_0 \f$. 
 *
 * As  this is multivariate normal, the the log-likelihood is then
 * \f[ -2 \log L = (X - E(X) )^T V^{-1} (X-E(X) ) + N\log(2\pi \det V) \f]
 *
 */
// Consider functionalizing this more so that most code can be reused in bm, bm_full, etc
// This should be able to store septime information, not call it each time the paramter values are changed!!
// i.e. septime should be passed to the algorithm, 
double ou_likelihood(tree * mytree)
{
	double sigma2 = mytree->pars[0];
	double alpha = mytree->pars[1];
	double theta = mytree->pars[2];

//	double root = mytree->trait[0];
	double root = theta;
	int n = (mytree->n+1)/2;
	int i, j;
	double s, t = 0;

	double * diff_X = (double *) malloc(n*sizeof(double));
	double * V = (double *) malloc(n*n*sizeof(double));
	int * tips = (int *) malloc(n * sizeof(int) );
	double * traits = (double *) malloc(n * sizeof(double) );

	gsl_matrix * V_inverse = gsl_matrix_alloc(n,n);
    gsl_permutation * p = gsl_permutation_alloc (n);
	gsl_matrix * ANS = gsl_matrix_alloc(1,n);
	gsl_matrix * ANS2 = gsl_matrix_alloc(1,1);

	double V_det, Xt_Vi_X;
	int signum;


	/* initialize just tips*/
	 j = 0;
	for(i=0; i < mytree->n; i++){
		if(mytree->tip[i])
		{ 
			tips[j] = i;
			traits[j] = mytree->trait[i];
			j++;
		}
	}

	/* calc total time in tree */
	i=0;
	while(!mytree->tip[i])
	{
		t += mytree->time[i];
		i = mytree->left[i];
	}
	t += mytree->time[i];


	/** @f$ E(X_t | X_0 ) = X_0 e^{-\alpha t} + \theta (1- e^{-\alpha t}) @f$ */
	double mean = (root-theta)*exp(-alpha*t) + theta;

	/* Compute E_x */
	for(i=0;i<n;i++)
	{
		diff_X[i] = traits[i] - mean;
	}

	/* Compute V: @f$ V_{ij} = \frac{\sigma^2}{2\alpha} (1 - e^{-2 \alpha s}) e^{-2\alpha (t-s)} @f$ */
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			s = t - sep_time(tips[i],tips[j],mytree)/2;
			V[n*i+j] = sigma2/(2*alpha)*(1-exp(-2*alpha*s) )* exp(-2*alpha*(t-s) );
		}
	}
	gsl_matrix_view V_view = gsl_matrix_view_array(V, n, n);
	gsl_matrix_view DIFF = gsl_matrix_view_array(diff_X, n, 1);
	gsl_linalg_LU_decomp (&V_view.matrix, p, &signum);
	gsl_linalg_LU_invert(&V_view.matrix, p, V_inverse);
	V_det = gsl_linalg_LU_det(&V_view.matrix,signum);


	/* @f$ -2 \log L = (X - E(X) )^T V^{-1} (X-E(X) ) + N\log(2\pi) + \log(\det V) @f$ */
	// Consider using appropriate blas optimized multiplication, not general matrix-matrix method!!
	gsl_blas_dgemm (CblasTrans, CblasNoTrans,
				   1.0, &DIFF.matrix, V_inverse,
				   0.0, ANS);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
				   1.0, ANS, &DIFF.matrix,
				   0.0, ANS2);
	Xt_Vi_X = gsl_matrix_get(ANS2, 0, 0);


	free(diff_X); 
	free(V);
	free(tips);
	free(traits);
	gsl_matrix_free(ANS);
	gsl_matrix_free(ANS2);
	gsl_matrix_free(V_inverse);
	gsl_permutation_free(p);

	return -Xt_Vi_X/2. -  n*log(2*M_PI)/2. - log(V_det)/2.;
}









/**
 *  Returns the log likelihood of brownian motion model (given tips, sigma2) 
 *
 *  The likelihood function under a linear model is a multivariate normal.  
 *  In the case of Brownian Motion, it is 
 *
 * */
double bm_likelihood(tree * mytree)
{
	double sigma2 = mytree->pars[0];

	double root = mytree->trait[0];
	int n = (mytree->n+1)/2;
	int i, j;
	double s, t = 0;

	double * diff_X = (double *) malloc(n*sizeof(double));
	double * V = (double *) malloc(n*n*sizeof(double));
	int * tips = (int *) malloc(n * sizeof(int) );
	double * traits = (double *) malloc(n * sizeof(double) );

	gsl_matrix * V_inverse = gsl_matrix_alloc(n,n);
    gsl_permutation * p = gsl_permutation_alloc (n);
	gsl_matrix * ANS = gsl_matrix_alloc(1,n);
	gsl_matrix * ANS2 = gsl_matrix_alloc(1,1);

	double V_det, Xt_Vi_X;
	int signum;


	/* initialize just tips*/
	 j = 0;
	for(i=0; i < mytree->n; i++){
		if(mytree->tip[i])
		{ 
			tips[j] = i;
			traits[j] = mytree->trait[i];
			j++;
		}
	}

	/* calc total time in tree */
	i=0;
	while(!mytree->tip[i])
	{
		t += mytree->time[i];
		i = mytree->left[i];
	}
	t += mytree->time[i];


	/* Compute E_x */
	for(i=0;i<n;i++)
	{
		diff_X[i] = traits[i] - root;
	}

	/* Compute V: @f$ V_{ij} =  @f$ */
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			s = t - sep_time(tips[i],tips[j],mytree)/2;
			V[n*i+j] = sigma2*s;
		}
	}
	gsl_matrix_view V_view = gsl_matrix_view_array(V, n, n);
	gsl_matrix_view DIFF = gsl_matrix_view_array(diff_X, n, 1);
	gsl_linalg_LU_decomp (&V_view.matrix, p, &signum);
	gsl_linalg_LU_invert(&V_view.matrix, p, V_inverse);
	V_det = gsl_linalg_LU_det(&V_view.matrix,signum);


	/* @f$ -2 \log L = (X - E(X) )^T V^{-1} (X-E(X) ) + N\log(2\pi \det V) @f$ */
	// Consider using appropriate blas optimized multiplication, not general matrix-matrix method!!
	gsl_blas_dgemm (CblasTrans, CblasNoTrans,
				   1.0, &DIFF.matrix, V_inverse,
				   0.0, ANS);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
				   1.0, ANS, &DIFF.matrix,
				   0.0, ANS2);
	Xt_Vi_X = gsl_matrix_get(ANS2, 0, 0);


	free(diff_X); 
	free(V);
	free(tips);
	free(traits);
	gsl_matrix_free(ANS);
	gsl_matrix_free(ANS2);
	gsl_matrix_free(V_inverse);
	gsl_permutation_free(p);

	return -Xt_Vi_X/2. -  n*log(2*M_PI)/2. - log(V_det)/2.;
}




/** Likelihood under BM of an ancestral state configuration (all internal as well as tip nodes) */
double bm_ancestral_likelihood(tree * mytree)
{
	return 1;
}




