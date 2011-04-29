/**
* @file likelihood.c
* @brief calculate the likelihood for general multitype OU processes on a phylogenetic tree
* @author Carl Boettiger <cboettig@gmail.com>
* @version 0.1
* @date 2011-04-22
*
* @section LICENSE 
*
* Copyright (C) 
* 2011 - Carl Boettiger
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
* 
*/

#include "likelihood.h"

/** get the last common ancestor of two nodes
 *  This version isn't particularly efficient, but the 
 *  calculation can in principle be done only once for 
 *  a given tree and stored as a matrix that is passed
 *  to the likelihood function.
 */
int 
get_lca(int i, int j, int n_nodes, const int * ancestor, 
            const double * branch_length, double * sep)
{
	int * ancestor_list = (int *) malloc(n_nodes * sizeof(int));
	int k = 0, s = 0;
	double * time_to_ancestor = (double *) calloc(n_nodes, sizeof(double));

	for(k=0; k<n_nodes; k++)
	{
		ancestor_list[k]=0;
	}
	k = 0;
	while(1){
		ancestor_list[k] = i;
		time_to_ancestor[k+1] = time_to_ancestor[k] + branch_length[i];
		if(i==0) break;
		i = ancestor[i];
		k++;
	}
	while(1){
		for(k=0; k<n_nodes; k++){
			if(j == ancestor_list[k]){ 
				s = j;
				j = 0;
				break; 
			}
		}
		if(j==0) break;
		j = ancestor[j];
	}
	*sep = time_to_ancestor[k];
	free(time_to_ancestor);
	free(ancestor_list);
	return s;	
}

/** Simple function to compute the age of a node 
 *  Called by calc_mean and calc_var */
double 
node_age(int i, const int * ancestor, const double * branch_length)
{
	double time=0;
	while(ancestor[i] >= 0 )
	{
		time += branch_length[i];
		i = ancestor[i];
	}
	return time; 
}

/**
* @brief Calculate the log normal likelihood given the vector
* of mean square differences and the variance matrix 
* @param n dimension 
* @param X_EX mean square differences, vector length n
* @param V variance matrix, n by n matrix (as a 1d array)
*
* @return Log likelihood
*/
double log_normal_lik(int n, double * X_EX, double * V)
{
	gsl_matrix * V_inverse = gsl_matrix_alloc(n,n);
    gsl_permutation * p = gsl_permutation_alloc (n);
	gsl_matrix * ANS = gsl_matrix_alloc(1,n);
	gsl_matrix * ANS2 = gsl_matrix_alloc(1,1);
	double V_det, Xt_Vi_X;
	int signum;

	gsl_matrix_view V_view = gsl_matrix_view_array(V, n, n);
	gsl_matrix_view DIFF = gsl_matrix_view_array(X_EX, n, 1);
	gsl_linalg_LU_decomp (&V_view.matrix, p, &signum);
	gsl_linalg_LU_invert(&V_view.matrix, p, V_inverse);
	V_det = gsl_linalg_LU_det(&V_view.matrix,signum);

	/** @f$ -2 \log L = (X - E(X) )^T V^{-1} (X-E(X) ) + N\log(2\pi) + \log(\det V) @f$ */
	// Consider using appropriate blas optimized multiplication, not general matrix-matrix method!!
	gsl_blas_dgemm (CblasTrans, CblasNoTrans,
				   1.0, &DIFF.matrix, V_inverse,
				   0.0, ANS);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
				   1.0, ANS, &DIFF.matrix,
				   0.0, ANS2);
	Xt_Vi_X = gsl_matrix_get(ANS2, 0, 0);

	gsl_matrix_free(ANS);
	gsl_matrix_free(ANS2);
	gsl_matrix_free(V_inverse);
	gsl_permutation_free(p);

	return -Xt_Vi_X/2. -  n*log(2*M_PI)/2. - log(V_det)/2.;
}


/** Calculate mean for likelihood computation, 
 */
double calc_mean(
	int i,
	double Xo, ///< root state
	const double * alpha, ///< value of alpha in each regime, length n_regimes
	const double * theta, ///< value of theta in each regime
	const int * regimes, ///< specification of the regimes (paintings), length n_nodes
	const int * ancestor, ///< ancestor of the node, length n_nodes
	const double * branch_length, ///< branch length ancestral to the node, length n_nodes
	double * output_gamma
	)
{
	int ri;
	/* Find the time until the root, which must have ancestor -1 */
	double time = node_age(i, ancestor, branch_length);
	double prev_time;

	/* Compute E(X_i) */
	long double gamma=0, omega=0;
	while( ancestor[i] >= 0 )
	{
		prev_time = time - branch_length[i]; 
		ri = regimes[i];
		gamma += alpha[ri] * branch_length[i];
		omega += theta[ri] * ( exp(alpha[ri]*time) - exp(alpha[ri] * prev_time) );
		i = ancestor[i];
		time = prev_time;
	}
	*output_gamma = gamma;

	return exp(-gamma)*(Xo + omega);
}



/**
 * Calculate the variance matrix for the multivariate normal likelihood
 * Function would be faster if it could reuse the gamma values calculated for the means.  
 * Function would be much faster if it stored gamma
 */
double calc_var(
	int i, int j, ///< nodes being compared
	int lca, ///< last common ancestor, pass to avoid recalculating
	const double * alpha, ///< value of alpha in each regime, length n_regimes
	const double * sigma, ///< value of sigma in each regime
	const int * regimes, ///< specification of the regimes (paintings), length n_nodes
	const int * ancestor, ///< ancestor of the node, length n_nodes
	const double * branch_length, ///< branch length ancestral to the node, length n_nodes
	const double * gamma_vec
	)
{
/** 
 * @f[ E(X_t) = \exp \left( - \sum \alpha_i \Delta t_i \right)
 * \left( X_{t_i} + \sum \theta_i \left(
 * e^{\alpha_i t_i}-e^{\alpha_i t_{i-1} } \right) \right) @f] */
	double gamma_i = gamma_vec[i], gamma_j = gamma_vec[j];

	double time = node_age(lca, ancestor, branch_length); 
	double prev_time;
	int ri;
	long double omega=0;

	i = lca;
	while( ancestor[i] >= 0) 
	{
		ri = regimes[i];
		prev_time = time - branch_length[i];
		omega += gsl_pow_2( sigma[ri] ) / (2 * alpha[ri] ) * (
			exp( 2 * alpha[ri] * time ) - exp( 2 * alpha[ri] * prev_time ) );
		time = prev_time;
		i = ancestor[i];
	}
	return exp(-gamma_i - gamma_j)*omega;
}


/**
* @brief Create a vector of length n tips with the identies of the tip nodes
*
* @param n_nodes Total number of nodes (internal and tips)
* @param ancestor List of ancestors
*
* @return an integer array of the tip nodes
*/
int * alloc_tips(int n_nodes, const int * ancestor){
	int n_tips = (n_nodes+1)/2;
	int * child = (int *) calloc(n_nodes,sizeof(int) );
	int * tips = (int *) calloc(n_tips,sizeof(int) );
	int i, j, k=0, empty;
	for(i=0; i < n_nodes; i++){ 
		empty = 1;
		for(j=0; j < n_nodes; j++){
			if(i == ancestor[j] ){
				if(empty){ 
					child[i] = j;
					empty = 0;
				}
			}
		}
		if (child[i] == 0){
			tips[k] = i;
			k++;
		}
	}
	free(child);
	return tips;
}




/**
* @brief Calculate the likelihood of a multitype OU model (callable from R) 
*
* @param Xo  Root value
* @param alpha[] selective strength in each regime 
* @param theta[] optimum in each regime
* @param sigma[] diversification rate in each regime
* @param regimes[] specify the regime each node belongs to
* @param ancestor[] specify the ancestor of each node (thus the topology)
* @param branch_length[] length of the branch below each node
* @param traits[] trait value observed at each node (tips only)
* @param n_nodes total number of nodes
* @param lca_matrix[] a n_nodes^2 matrix of least common ancestor for each pair,
*  computed by the lca_calc function once per tree for efficiency.
* @param llik the likelihood returned by the function
*/
void calc_lik (const double *Xo, const double alpha[], const double theta[], 
	             const double sigma[], const int regimes[], const int ancestor[],
                 const double branch_length[], const double traits[], 
                 int *n_nodes, int lca_matrix[], double *llik)
{
	gsl_set_error_handler_off ();

	int i,j,ki, kj;
	int n_tips = (*n_nodes+1)/2;
	double *X_EX = (double *) malloc(n_tips * sizeof(double));
	double *V = (double *) malloc(n_tips * n_tips * sizeof(double));
	double *gamma_vec = (double *) calloc(*n_nodes, sizeof(double));
	double mean, gamma_i;
	int lca;

	int * tips = alloc_tips(*n_nodes, ancestor);

    /* Calculate the mean square differences */
	for(i = 0; i < n_tips; i++){
		ki = tips[i];
		mean = calc_mean(ki, *Xo, alpha, theta, regimes, ancestor,
                         branch_length, &gamma_i);
		X_EX[i] = traits[ki] - mean;
		gamma_vec[ki] = gamma_i;
	}
    /* Calculate the varainces */
	for(i=0; i < n_tips; i++){
		ki = tips[i];
		for(j=0; j< n_tips; j++){
			kj = tips[j];
			lca = lca_matrix[ki * *n_nodes+kj];
			V[n_tips*i+j] = calc_var(ki,kj,lca, alpha, sigma, regimes,
                                     ancestor, branch_length, gamma_vec);
		}
	}


	*llik = log_normal_lik(n_tips, X_EX, V);
	free(gamma_vec);
	free(X_EX); 
	free(V);
	free(tips);
}


/**
* @brief Simulate by drawing random numbers
* from the appropriate multivariate normal 
*
* @param rng a gsl random number gerneator
* @param mytree the tree and parameters
*/
void simulate (const gsl_rng * rng, tree * mytree)
{
  int i,j,ki, kj;
  int n_tips = (mytree->n_nodes+1)/2;
  gsl_vector * EX = gsl_vector_alloc(n_tips);
  gsl_matrix * V = gsl_matrix_alloc(n_tips,n_tips);
  gsl_vector * simdata = gsl_vector_alloc(n_tips);
  double * gamma_vec = (double *) calloc(mytree->n_nodes,sizeof(double));
  int * tips = alloc_tips(mytree->n_nodes, mytree->ancestor);

  double gamma_i;
  int lca;

  for(i = 0; i < n_tips; i++){
    ki = tips[i];
    gsl_vector_set( EX,
                    i,
                    calc_mean(ki, *(mytree->Xo), mytree->alpha, 
                              mytree->theta, mytree->regimes, 
                              mytree->ancestor, mytree->branch_length, 
                              &gamma_i));
    gamma_vec[ki] = gamma_i;
  }
  for(i=0; i < n_tips; i++){
    ki = tips[i];
    for(j=0; j< n_tips; j++){
      kj = tips[j];
      lca = mytree->lca_matrix[ki*mytree->n_nodes+kj];
      gsl_matrix_set(	V, i, j, 
                      calc_var(ki,kj,lca, mytree->alpha, 
                               mytree->sigma, mytree->regimes, 
                               mytree->ancestor, mytree->branch_length, 
                               gamma_vec));
    }
  }

  mvn(rng, EX, V, simdata);

  for(i=0; i< n_tips; i++){
      ki = tips[i];
      mytree->traits[ki] = gsl_vector_get(simdata,i);
  }

  gsl_vector_free(EX);
  gsl_matrix_free(V);
  gsl_vector_free(simdata);
  free(gamma_vec);
  free(tips);
  }



