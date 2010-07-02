#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

int get_lca(int i, int j, int n_nodes, const int * ancestor)
	{
	int * ancestor_list = (int *) malloc(n_nodes * sizeof(int));
	int k = 0, s = 0;
	for(k=0; k<n_nodes; k++)
	{
		ancestor_list[k]=0;
	}
	k = 0;
	while(1){
		ancestor_list[k] = i;
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
	free(ancestor_list);
	return s;	
}




double node_age(int i, const int * ancestor, const double * branch_length){
	double time=0;
	while(ancestor[i] >= 0 )
	{	
		time += branch_length[i];
		i = ancestor[i];
	}
	return time; 
}

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

	/* @f$ -2 \log L = (X - E(X) )^T V^{-1} (X-E(X) ) + N\log(2\pi) + \log(\det V) @f$ */
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


/**
 * @f[ E(X_t) = \exp \left( - \sum \alpha_i \Delta t_i \right) \left( X_) + \sum \theta_i \left( e^{\alpha_i t_i}-e^{\alpha_i t_{i-1} } \right) \right)
 */
double calc_mean(
	int i,
	double Xo, ///< root state
	const double * alpha, ///< value of alpha in each regime, length n_regimes
	const double * theta, ///< value of theta in each regime
	const double * sigma, ///< value of sigma in each regime
	const int * regimes, ///< specification of the regimes (paintings), length n_nodes
	const int * ancestor, ///< ancestor of the node, length n_nodes
	const double * branch_length ///< branch length ancestral to the node, length n_nodes
	)
{
	int ri;
	/* Find the time until the root, which must have ancestor -1 */
	double time = node_age(i, ancestor, branch_length);
	double prev_time;

	/* Compute E(X_i) */
	double gamma=0, omega=0;
	while( ancestor[i] >= 0 )
	{
		prev_time = time - branch_length[i]; 
		ri = regimes[i];
		gamma += alpha[ri] * branch_length[i];
		omega += theta[ri] * ( exp(alpha[ri]*time) - exp(alpha[ri] * prev_time) );
		i = ancestor[i];
		time = prev_time;
	}
	return exp(-gamma)*(Xo + omega);
}



double calc_gamma(int i, const int * ancestor, const double * branch_length, const int * regimes, const double * alpha){
	double gamma = 0;
	while( ancestor[i] >= 0 )
	{
		gamma += alpha[regimes[i]] * branch_length[i];
		i = ancestor[i];
	}
	return gamma;
}

/**
 * @f[ E(X_t) = \exp \left( - \sum \alpha_i \Delta t_i \right) \left( X_) + \sum \theta_i \left( e^{\alpha_i t_i}-e^{\alpha_i t_{i-1} } \right) \right)
 */
double calc_var(
	int i, int j, ///< nodes being compared
	int lca, ///< last common ancestor
	double Xo, ///< root state
	const double * alpha, ///< value of alpha in each regime, length n_regimes
	const double * theta, ///< value of theta in each regime
	const double * sigma, ///< value of sigma in each regime
	const int * regimes, ///< specification of the regimes (paintings), length n_nodes
	const int * ancestor, ///< ancestor of the node, length n_nodes
	const double * branch_length ///< branch length ancestral to the node, length n_nodes
	)
{
	double gamma_i = calc_gamma(i, ancestor, branch_length, regimes, alpha);
	double gamma_j = calc_gamma(j, ancestor, branch_length, regimes, alpha);
	double time = node_age(lca, ancestor, branch_length); 
	double prev_time;
	int ri;
	double omega=0;

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


int * alloc_tips(int n_nodes, const int * ancestor){
	int n_tips = (n_nodes+1)/2;
	int * child = (int *) malloc(n_nodes*sizeof(int) );
	int * tips = (int *) malloc(n_tips*sizeof(int) );
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


double gen_lik(	double Xo, ///< root state
	const double * alpha, ///< value of alpha in each regime, length n_regimes
	const double * theta, ///< value of theta in each regime
	const double * sigma, ///< value of sigma in each regime
	const int * regimes, ///< specification of the regimes (paintings), length n_nodes
	const int * ancestor, ///< ancestor of the node, length n_nodes
	const double * branch_length, ///< branch length ancestral to the node, length n_nodes
	const double * traits, ///< traits
	int n_nodes
	)
{
	int i,j,ki, kj;
	int n_tips = (n_nodes+1)/2;
	double * X_EX = (double *) malloc(n_tips*sizeof(double));
	double * V = (double *) malloc(n_tips*n_tips*sizeof(double));

	double llik;
	int lca;

	int * tips = alloc_tips(n_nodes, ancestor);
	
	for(i = 0; i < n_tips; i++){
		ki = tips[i];
		X_EX[i] = traits[ki] - calc_mean(ki, Xo, alpha, theta, sigma, regimes, ancestor, branch_length);
	}
	for(i=0; i < n_tips; i++){
		ki = tips[i];
		for(j=0; j< n_tips; j++){
			kj = tips[j];
			lca = get_lca(ki,kj, n_nodes, ancestor);
			V[n_tips*i+j] = calc_var(ki,kj,lca, Xo, alpha, theta, sigma, regimes, ancestor, branch_length);
//			printf("%.2lf ", V[n_tips*i+j]);
		}
//		printf("\n");
	}
//	printf("\n");

	llik = log_normal_lik(n_tips, X_EX, V);
	free(X_EX); 
	free(V);
	free(tips);
	return llik;
}

int main(void)
{
	/*
	double branch_length[] =  {0, 1, 1, 3, 2, 1, 1}; 
	const int ancestor[] = {-1, 2, 0, 0, 2, 1, 1};
	const double traits[] = {0, 0, 0, -2, 1, 2, 2};	
	const int regimes[] =	{0, 0, 0, 0, 0, 0, 0};
	int n_nodes = 7;
*/
	int n_nodes = 45;
	double branch_length[] = { 0, 12./38, 20./38,  2./38,  2./38,  4./38,  
							 8./38,  5./38,  5./38,  5./38,  5./38, 10./38,  
							 9./38,  4./38,  8./38,  2./38, 20/38,  2./38,  
							 4./38,  2./38, 1./38,  2./38, 26./38, 4./38, 
							 2./38, 2./38,  2./38,  2./38, 15./38, 10./38, 
							 10./38, 10./38, 10./38, 16./38, 12./38,  4./38,  
							 2./38,  2./38, 10./38,  8./38,  2./38,  1./38,  
							 1./38,  2./38,  2./38};

	int ancestor[] = {-1,  0,  1,  2,  3,  2,  0,  6,  7,  8,  9,
						8,  7, 12, 13,14,  6, 16, 17, 18, 19, 18, 1,
						3,  4,  4,  5,  5,  9, 10, 10, 11, 11, 12, 
						13, 14, 15, 15, 16, 17, 19, 20, 20, 21, 21};

	double traits[] = { 2.900000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
						  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
						  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.602690, 2.660260,
						  2.660260, 2.653242, 2.674149, 2.701361, 3.161247, 3.299534, 3.328627, 3.353407,
						  3.360375, 3.049273, 2.906901, 2.980619, 2.933857, 2.975530, 3.104587, 3.346389,
						  2.928524, 2.939162, 2.990720, 3.058707, 3.068053};
	int regimes[45] = {0};
	int LP_regimes[] =  {1, 1, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
						2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1}; 


	double Xo = 2.9;
	double alpha[1] = {0.1921554};
	double theta[1] = {2.953806};
	double sigma[1] = {sqrt(0.04836469)};

	printf("-llik = %lf\n", gen_lik(Xo, alpha, theta, sigma, regimes, ancestor, branch_length, traits, n_nodes) );
/*
	int i = 3; int j = 3;

	double t = node_age(i, ancestor, branch_length);
	double mean = (Xo-theta[0])*exp(-alpha[0]*t) + theta[0];
	double Ex = calc_mean(i, Xo, alpha, theta, sigma, regimes, ancestor, branch_length);
	printf("Mean: %lf, %lf\n", Ex, mean);

	// double s = 2;
	//double covar = gsl_pow_2(sigma[0]) / (2*alpha[0]) * ( 1 - exp(-2*alpha[0] * s) )*exp(-2*alpha[0]*(t-s) );
	
	int lca = get_lca(i,j, n_nodes, ancestor);
	double Vx = calc_var(i, j, lca, Xo, alpha, theta, sigma, regimes, ancestor, branch_length);
	printf("Var(%d,%d): %lf\n",i,j, Vx);
*/

	return 0;
}
