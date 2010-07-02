#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
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

	/* Find the current time of each node (redundant if both are tips) */
	double time_i = node_age(i, ancestor, branch_length);
	double time_j = node_age(j, ancestor, branch_length);

	double gamma_i = calc_gamma(i, ancestor, branch_length, regimes, alpha);
	double gamma_j = calc_gamma(j, ancestor, branch_length, regimes, alpha);

	double prev_time_i = time_i, prev_time_j = time_j;
	int ri, rj;
	double omega=0;

	/* Calculate contribution to the variance/covariance since divergence 
	while( i != lca && j != lca)
	{
		ri = regimes[i];
		rj = regimes[j];
		prev_time_i = time_i - branch_length[i];
		prev_time_j = time_j - branch_length[j];
		omega += sigma[ri] * sigma[rj] / ( alpha[ri]+alpha[rj] ) * ( 
			exp(  (alpha[ri] + alpha[rj]) * GSL_MIN( time_i, time_j) ) - 
			exp(  (alpha[ri] + alpha[rj]) * GSL_MAX( prev_time_i, prev_time_j))  );

		if(prev_time_i >= prev_time_j) {
			i = ancestor[i];
			time_i = prev_time_i;
		} else if (prev_time_i < prev_time_j) {
			time_j = prev_time_j;
			j = ancestor[j];
		}
	}
	*/

	/* Calculate contribution to the variance/covariance prior to divergence */
	i = lca;
	time_i = node_age(i, ancestor, branch_length); // should also equal prev_time_j
	while( ancestor[i] >= 0) // Note i = j = LCA
	{
printf("shared\n");

		ri = regimes[i];
		prev_time_i = time_i - branch_length[i];
		omega += gsl_pow_2( sigma[ri] ) / (2 * alpha[ri] ) * (
			exp( 2 * alpha[ri] * time_i ) - exp( 2 * alpha[ri] * prev_time_i ) );
		time_i = prev_time_i;
		i = ancestor[i];
	}
printf("e= %lf, %lf\n", exp(-gamma_i - gamma_j), gamma_i+gamma_j);
	return exp(-gamma_i - gamma_j)*omega;
}

int main(void)
{
	// test data
//	const double traits[] = {0, 0, 0, -2, 1, 2, 2};
	
	const double branch_length[7] =  {0, 1, 1, 3, 2, 1, 1}; 
	const int ancestor[7] = {-1, 2, 0, 0, 2, 1, 1};
	const int regimes[7] =	{0, 0, 0, 0, 0, 0, 0};

	double Xo = 0;
	double alpha[1] = {1};
	double theta[1] = {2};
	double sigma[1] = {1};

	int i = 4; int j = 6;

	double t = node_age(i, ancestor, branch_length);

	double mean = (Xo-theta[0])*exp(-alpha[0]*t) + theta[0];
	double gamma = calc_gamma(i, ancestor, branch_length, regimes, alpha);
	printf("%lf %lf \n", alpha[0]*t, gamma);

	double Ex = calc_mean(i, Xo, alpha, theta, sigma, regimes, ancestor, branch_length);
	printf("%lf, %lf\n", Ex, mean);

	double s = 1;
	double covar = gsl_pow_2(sigma[0]) / (2*alpha[0]) * ( 1 - exp(-2*alpha[0] * s) )*exp(-2*alpha[0]*(t-s) );

	int lca = get_lca(i,j, 7, ancestor);
	printf("lca = %d, t=%lf\n", lca, t);

	printf("%lf %lf\n", calc_var(i, j, lca, Xo, alpha, theta, sigma, regimes, ancestor, branch_length), covar);

	return 0;
}
