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

/** Simple function to compute the age of a node measured from root, forward 
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

  int is_singular; /* error checking */

  gsl_matrix_view V_view = gsl_matrix_view_array(V, n, n);
  gsl_matrix_view DIFF = gsl_matrix_view_array(X_EX, n, 1);
  is_singular = gsl_linalg_LU_decomp (&V_view.matrix, p, &signum);
  if(is_singular) return GSL_NEGINF;
  gsl_linalg_LU_invert(&V_view.matrix, p, V_inverse);
  V_det = gsl_linalg_LU_det(&V_view.matrix,signum);

  /**@f$ -2 log L = @f$  
   * @f$ (X - E(X) )^T V^{-1} (X-E(X) ) + N\log(2\pi) + \log(\det V) @f$ */
  /* Consider using appropriate blas optimized multiplication, 
   * not the general matrix-matrix method (for greater speed) */
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

/** Calculate the row of the gamma matrix for tip i
 * @param tips a vector giving the node number of each tip (data)
 * @param n_tips the number of tips (data)
 * @param alpha vector of alpha parameters for each regime (data)
 * @param regimes the regime of each branch (data)
 * @param branch_length vector of branch lengths of each node (data)
 * @param gamma_matrix the returned matrix (well, a pointer to it)
 * @details in gamma_matrix[i,j], i is the node-number of the tip
 * whose history we are tracking, j is any other node.  Only for j's
 * in the ancestry of i will gamma[i,j] have entries.  Those entries
 * are simply the sum of branch lengths times alpha values of that node */
void calc_gamma_matrix(const int * tips, const int n_tips, 
                       const double * alpha, const int * regimes, 
                       const int * ancestor, const double * branch_length,
                       gsl_matrix * gamma_matrix)
{
  int i, j, k, rj;
  double value;
  for(k = 0; k < n_tips; k++)
  {
    i = tips[k]; 
    j = i;
    value = 0;
    while( ancestor[j] >= 0)
    {
      rj = regimes[j];
      value += branch_length[j] * alpha[rj];
      j = ancestor[j];
      gsl_matrix_set(gamma_matrix, i, j, value);
    }
  }
}






 
double calc_mean(int i, double Xo,  const double * alpha, 
                 const double * theta, const int * regimes, 
                 const int * ancestor, const double * branch_length, 
                 const gsl_matrix * gamma_matrix)
{
  int tip = i; // remember the tip value
  int ri;
  long double omega=0;
  while( ancestor[i] >= 0 )
  {
    ri = regimes[i];
    omega += theta[ri] * (1 - exp(- alpha[ri] * branch_length[i])) *
             exp( - gsl_matrix_get(gamma_matrix, tip, i) );
    i = ancestor[i];
  }
  return Xo * exp( - gsl_matrix_get(gamma_matrix, tip, 0)) + omega;
}



/** Calculate the covariance between tips i & j
 * @f[ e^{-\gamma_{i,k}} e^{-\gamma_{j_k}} \frac{\sigma_k^2}{2 \alpha_k} 
 *     \left(e^{2\alpha_k T } - e^{2 \alpha_k (T - l)} \right) @f]
 *  where T is the age of the node and l the length of the branch under it
 */
double calc_var(
  int i, int j, ///< nodes being compared
  int lca, ///< last common ancestor, pass to avoid recalculating
  const double * alpha, ///< value of alpha in each regime, length n_regimes
  const double * sigma, ///< value of sigma in each regime
  const int * regimes, ///< specification of the regimes (paintings), length n_nodes
  const int * ancestor, ///< ancestor of the node, length n_nodes
  const double * branch_length, ///< branch length ancestral to the node, length n_nodes
  gsl_matrix * gamma_matrix 
  )
{
  int k = lca; //if i=j, k=i=j 
  int rk;
  long double omega=0;
	
  while( ancestor[k] >= 0) 
  {
    rk = regimes[k];
    omega += gsl_pow_2( sigma[rk] ) / (2 * alpha[rk] ) * 
      ( 1 - exp( - 2 * alpha[rk] * branch_length[k] ) ) *
      exp( - gsl_matrix_get(gamma_matrix, i, k) - gsl_matrix_get(gamma_matrix, j, k));
    k = ancestor[k];
  }
  return omega;
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
*
* @details 
*  Calculates the expected value and the convariance matrix for each tip.  
*
*/
void calc_lik (const double *Xo, const double alpha[], const double theta[], 
               const double sigma[], const int regimes[], const int ancestor[],
                 const double branch_length[], const double traits[], 
                 int *n_nodes, int lca_matrix[], double *llik)
{
  /* gsl_set_error_handler_off (); /* Comment out this line to assist debugging */

  /* Declare variables */
  int i, j, ki, kj;
  int n_tips = (*n_nodes+1)/2;
  double *X_EX = (double *) malloc(n_tips * sizeof(double));
  double *V = (double *) malloc(n_tips * n_tips * sizeof(double));
  gsl_matrix * gamma_matrix = gsl_matrix_calloc(*n_nodes,*n_nodes);
  double mean;
  int lca;
  int * tips = alloc_tips(*n_nodes, ancestor);


  /* Calculate the gamma matrix */
  calc_gamma_matrix(tips, n_tips, alpha, regimes, ancestor,
                    branch_length, gamma_matrix);


  /* Unit test -- tips have the same age */
//  for(i = 0; i < n_tips; i++)
//    printf("%lf\n", node_age(tips[i], ancestor, branch_length));

  /* Unit test -- gamma of root values */
 // for(i = 0; i < n_tips; i++)
//    printf("%lf\n", gsl_matrix_get(gamma_matrix, tips[i], 0));
//    printf("%lf\n", gsl_matrix_get(gamma_matrix, tips[i], ancestor[tips[i]]));

/* Unit test -- variance on a single tip with a middle node (tree = *-*-*)  */
/*
double salpha[] = {.1}; 
double ssigma[] = {2}; 
int sregimes[] = {0, 0, 0, 0}; 
int sancestor[] = {-1, 0, 1, 1}; 
double sbranch_length[] = {0, 5, 5, 5}; 
int s_tips[] = {2,3};
gsl_matrix * sgamma_matrix = gsl_matrix_calloc(4,4);
calc_gamma_matrix(s_tips, 2, salpha, sregimes, sancestor,
                    sbranch_length, sgamma_matrix);
printf("var: %lf\n",	calc_var(2, 3, 1, salpha, ssigma, sregimes, sancestor, sbranch_length, sgamma_matrix));
printf("analytic %lf\n", gsl_pow_2(ssigma[0])/(2*salpha[0]) * (1 - exp(-2*salpha[0]*5)) *exp(-2*salpha[0]*5) );

printf("var: %lf\n",	calc_var(2, 2, 2, salpha, ssigma, sregimes, sancestor, sbranch_length, sgamma_matrix));
printf("analytic %lf\n", gsl_pow_2(ssigma[0])/(2*salpha[0]) * (1 - exp(-2*salpha[0]*10)));
*/
/*
for(i=0;i<4;i++){
	printf("\n");
  for(j=0;j<4;j++){
		printf("%g\t ", gsl_matrix_get(sgamma_matrix, i, j));
	}
}
printf("\n\n");
*/


/*
for(i=0;i<*n_nodes;i++){
	printf("\n");
  for(j=0;j<15;j++){
		printf("%g\t ", gsl_matrix_get(gamma_matrix, i, j));
	}
}
printf("\n\n");
*/

  /* Calculate the mean square differences */
  for(i = 0; i < n_tips; i++){
    ki = tips[i];
    mean = calc_mean(ki, *Xo, alpha, theta, regimes, ancestor,
                     branch_length, gamma_matrix);
    X_EX[i] = traits[ki] - mean;
//printf("%lf\n", mean);
  }


  /* Calculate the variances */
  for(i=0; i < n_tips; i++){
    ki = tips[i];
    for(j=0; j < n_tips; j++){
      kj = tips[j];
      /* Identify which node is last common ancestor of the tips*/
      lca = lca_matrix[ki * *n_nodes + kj];
      /* get the covariance between all possible pairs of tips */
      V[n_tips*i+j] = calc_var(ki, kj, lca, alpha, sigma, regimes,
                               ancestor, branch_length, gamma_matrix);
//if(ki==kj) printf("%g, %d, %d\n", V[n_tips*i+j], ki, lca);
    }
  }

  *llik = log_normal_lik(n_tips, X_EX, V);
  gsl_matrix_free(gamma_matrix);
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
  /* Allocate memory */
  int i,j,ki, kj;
  int n_tips = (mytree->n_nodes+1)/2;
  gsl_vector * EX = gsl_vector_alloc(n_tips);
  gsl_matrix * V = gsl_matrix_alloc(n_tips,n_tips);
  gsl_vector * simdata = gsl_vector_alloc(n_tips);
  gsl_matrix * gamma_matrix = gsl_matrix_alloc(mytree->n_nodes,mytree->n_nodes);
  int * tips = alloc_tips(mytree->n_nodes, mytree->ancestor);
  int lca;

  /* Calculate the gamma matrix */
  calc_gamma_matrix(tips, n_tips, mytree->alpha, mytree->regimes,
                    mytree->ancestor, mytree->branch_length, 
                    gamma_matrix);

  /* Calculate means */
  for(i = 0; i < n_tips; i++){
    ki = tips[i];
    gsl_vector_set(EX, i,
                    calc_mean(ki, *(mytree->Xo), mytree->alpha, 
                              mytree->theta, mytree->regimes, 
                              mytree->ancestor, mytree->branch_length, 
                              gamma_matrix));
  }
  /* Calculate Variances */
  for(i=0; i < n_tips; i++){
    ki = tips[i];
    for(j=0; j< n_tips; j++){
      kj = tips[j];
      lca = mytree->lca_matrix[ki * mytree->n_nodes + kj];
      gsl_matrix_set( V, i, j, 
                      calc_var(ki,kj,lca, mytree->alpha, 
                               mytree->sigma, mytree->regimes, 
                               mytree->ancestor, mytree->branch_length, 
                               gamma_matrix));
    }
  }

  /* Calculate simulated data as multivariate normal random numbers */
  mvn(rng, EX, V, simdata);

  /* Write that data to the tip states */
  for(i=0; i< n_tips; i++){
      ki = tips[i];
      mytree->traits[ki] = gsl_vector_get(simdata,i);
  }

  /* Clean up */
  gsl_vector_free(EX);
  gsl_matrix_free(V);
  gsl_vector_free(simdata);
  gsl_matrix_free(gamma_matrix);
  free(tips);
}


void unit_tests (const double *Xo, const double alpha[], const double theta[], 
               const double sigma[], const int regimes[], const int ancestor[],
                 const double branch_length[], const double traits[], 
                 int *n_nodes, int lca_matrix[], double *llik)
{
  /* gsl_set_error_handler_off (); /* Comment out this line to assist debugging */

  /* Declare variables */
  int i, j, ki, kj;
  int n_tips = (*n_nodes+1)/2;
  double *X_EX = (double *) malloc(n_tips * sizeof(double));
  double *V = (double *) malloc(n_tips * n_tips * sizeof(double));
  gsl_matrix * gamma_matrix = gsl_matrix_calloc(*n_nodes,*n_nodes);
  double mean;
  int lca;
  int * tips = alloc_tips(*n_nodes, ancestor);


  /* Calculate the gamma matrix */
  calc_gamma_matrix(tips, n_tips, alpha, regimes, ancestor,
                    branch_length, gamma_matrix);


  /* Unit test -- tips have the same age */
  for(i = 0; i < n_tips; i++)
    printf("%lf\n", node_age(tips[i], ancestor, branch_length));

  /* Unit test -- gamma of root values */
  for(i = 0; i < n_tips; i++)
    printf("%lf\n", gsl_matrix_get(gamma_matrix, tips[i], 0));
    printf("%lf\n", gsl_matrix_get(gamma_matrix, tips[i], ancestor[tips[i]]));

/* Unit test -- variance on a single tip with a middle node (tree = *-*-*)  */
double salpha[] = {.1}; 
double ssigma[] = {2}; 
int sregimes[] = {0, 0, 0, 0}; 
int sancestor[] = {-1, 0, 1, 1}; 
double sbranch_length[] = {0, 5, 5, 5}; 
int s_tips[] = {2,3};
gsl_matrix * sgamma_matrix = gsl_matrix_calloc(4,4);
calc_gamma_matrix(s_tips, 2, salpha, sregimes, sancestor,
                    sbranch_length, sgamma_matrix);
printf("var: %lf\n",	calc_var(2, 3, 1, salpha, ssigma, sregimes, sancestor, sbranch_length, sgamma_matrix));
printf("analytic %lf\n", gsl_pow_2(ssigma[0])/(2*salpha[0]) * (1 - exp(-2*salpha[0]*5)) *exp(-2*salpha[0]*5) );

printf("var: %lf\n",	calc_var(2, 2, 2, salpha, ssigma, sregimes, sancestor, sbranch_length, sgamma_matrix));
printf("analytic %lf\n", gsl_pow_2(ssigma[0])/(2*salpha[0]) * (1 - exp(-2*salpha[0]*10)));

for(i=0;i<4;i++){
	printf("\n");
  for(j=0;j<4;j++){
		printf("%g\t ", gsl_matrix_get(sgamma_matrix, i, j));
	}
}
printf("\n\n");


  /* Calculate the mean square differences */
  for(i = 0; i < n_tips; i++){
    ki = tips[i];
    mean = calc_mean(ki, *Xo, alpha, theta, regimes, ancestor,
                     branch_length, gamma_matrix);
    X_EX[i] = traits[ki] - mean;
    printf("%lf\n", mean);
  }


  /* Calculate the variances */
  for(i=0; i < n_tips; i++){
    ki = tips[i];
    for(j=0; j < n_tips; j++){
      kj = tips[j];
      /* Identify which node is last common ancestor of the tips*/
      lca = lca_matrix[ki * *n_nodes + kj];
      /* get the covariance between all possible pairs of tips */
      V[n_tips*i+j] = calc_var(ki, kj, lca, alpha, sigma, regimes,
                               ancestor, branch_length, gamma_matrix);
      if(ki==kj) printf("%g, %d, %d\n", V[n_tips*i+j], ki, lca);
    }
  }

  *llik = log_normal_lik(n_tips, X_EX, V);
  gsl_matrix_free(gamma_matrix);
  free(X_EX); 
  free(V);
  free(tips);
}


