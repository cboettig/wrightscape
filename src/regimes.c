#include "tree.h"

/** The stationary distribution of the OU process */
double stationary(int s, double x, tree * mytree)
{
	double sigma2= mytree->pars[3*s];
	double alpha = mytree->pars[3*s+1];
	double theta = mytree->pars[3*s+2];
	return exp(-gsl_pow_2(x-theta)*alpha/(sigma2))/sqrt(4*M_PI*sigma2/alpha); 
}

/** The time dependent solution of the OU process under specified regime */
double regimes_ou(int si, double xi, int sj, double xj, double t, tree * mytree)
{
	double sigma2= mytree->pars[3*si];
	double alpha = mytree->pars[3*si+1];
	double theta = mytree->pars[3*si+2];
	double u = (1-exp(-alpha*t) ) * theta + xi*exp(-alpha*t);
	double v = (1-exp(-2*alpha*t) )* sigma2/(2*alpha) ;
	double ans = exp(-gsl_pow_2( xj-u ) / (2*v) ) / sqrt(2*M_PI*v);
	return ans;
}


/** 
 * Create a stochastic transition matrix from a vector of the rates:
 * i.e. the vector \f$ \lambda_{12}, \lambda_{13}, \lambda_{21}, \lambda_{23}, \lambda{31}, \lambda{32} \f$
 * becomes the matrix:
 * \f[ \begin{pmatrix}
 *		-(\lambda_{12}+\lambda_{13}) & \lambda{12}					 & \lambda{13} \\
 *		\lambda_{21}				 & -(\lambda_{21} +\lambda_{23} ) & \lambda{23} \\
 *		\lambda_{31}				 & \lambda_{32}					 & -(\lambda_{31}+\lambda_{32} ) \\
 *	\end{pmatrix} \f]
 *
 * */
gsl_matrix * rates_to_Q(tree * mytree)
{
	size_t n_regimes = mytree->nstates;
	gsl_matrix * Q = gsl_matrix_alloc(n_regimes, n_regimes);
	int i,j; 
	for(i=0;i<n_regimes; i++){
		double sum = 0;
		for(j=0; j<n_regimes-1; j++){
			sum += mytree->pars[i*(n_regimes-1)+j+3*n_regimes];
			gsl_matrix_set(
				Q, 
				i, 
				(1+j+i)%n_regimes, 
				mytree->pars[i*(n_regimes-1)+j+3*n_regimes]
			);
		}
		gsl_matrix_set(Q, i, i, -sum); 
	}	
	return Q;
}


/** jump from si, xi to sj, xj in time t */
double regimes(int si, double xi, int sj, double xj, double t, tree * mytree)
{
	size_t n_regimes = mytree->nstates;
	gsl_matrix * Q = rates_to_Q(mytree);
	gsl_matrix * expQ = gsl_matrix_alloc(n_regimes, n_regimes);
	gsl_matrix_scale(Q, t);
	gsl_linalg_exponential_ss(Q, expQ, .01);

	double Pi = stationary(sj, xj, mytree);

	double P_in_i = exp( gsl_matrix_get(Q, si, si)*t );
	double P_jumped = gsl_matrix_get(expQ, si, sj)*Pi;
	double P_no_jump = P_in_i*Pi;
	double P_stay = P_in_i * regimes_ou(si, xi, sj, xj, t, mytree); 


	gsl_matrix_free(Q);
	gsl_matrix_free(expQ);

	return P_jumped + (si==sj) * (- P_no_jump +P_stay );
}

/** 
 * All matrices are stored in one big vector.  Organized as follows:
 *
 * Each node gets a continuous block in order of their numbering.  
 *
 * Within their block, tips have a block of length grid-size for the 
 * transition probabilities coming from each state, in order.  Hence
 * the size of a block for a tip is gridsize*n_states.
 * Hence the matrix for a tip is located at node_begins + s*gridsize
 *
 * Internal nodes are size gridsize^2 * n_states^2
 * a matrix for internal node from state si to sj is found at:
 * node_begins+s1*n_states*gsl_pow_2(gridsize) + sj*gsl_pow_2(gridsize)
 *
 * */
int * indices_alloc_init(tree * mytree)
{
	size_t n = mytree->n, n_states = mytree->nstates, gridsize = mytree->gridsize;
	int * node_begins = (int *) malloc(n*sizeof(int));
	int i, node_block, current_position = 0;
	for(i=0;i< n; i++){
		if(mytree->tip[i] || i == 0){ 
			node_block =  gridsize*n_states;
		} else {
			node_block = gsl_pow_2(gridsize)*gsl_pow_2(n_states);
		}
		node_begins[i] = current_position;
		current_position += node_block;
	}
	return node_begins;
}

/**
 * initialize the matrices for all branches under all regimes.  
 */
void regimes_matrices_init(gsl_vector * regimes_matrices, tree * mytree){
	int node=0, si, sj, i, j, index;
	size_t n = mytree->n, n_states = mytree->nstates, gridsize = mytree->gridsize;
	int * node_begins = indices_alloc_init(mytree);


	/* Set up the grid */
	double delta_grid = (double) (ULIM-LLIM)/(gridsize-1);
	double * x = (double *) malloc(gridsize*sizeof(double) );
	for(i=0;i<gridsize;i++){
		x[i] = (double) delta_grid*i + LLIM;
	}

	/* root */
	int root_flag = 1;
	for(si=0;si<n_states; si++){
		index = node_begins[node] + si * gridsize;
		for(i=0;i<gridsize;i++){ 
			gsl_vector_set(regimes_matrices, index, 0);
			if(i*delta_grid+LLIM > mytree->trait[0] && root_flag){
				gsl_vector_set(regimes_matrices, index, 1);
				root_flag = 0;
			}
			index++;
		}
	}
	for(node=1;node<n;node++){
		if(mytree->tip[node]){ 
			/* tips */
			for(sj=0;sj<n_states; sj++){
				index = node_begins[node] + sj * gridsize;
				for(j=0;j<gridsize;j++){
					gsl_vector_set(regimes_matrices, index, regimes(mytree->state[node], mytree->trait[node],  sj, x[j], mytree->time[node], mytree) );
					index++;
				}
			}
		} else {
			/* internal nodes */
			for(si=0;si<n_states; si++){
				for(sj=0;sj<n_states; sj++){
					index = node_begins[node] + si * n_states * gsl_pow_2(gridsize) + sj * gsl_pow_2(gridsize);
					for(i=0;i<gridsize;i++){
						for(j=0;j<gridsize;j++){
							gsl_vector_set(regimes_matrices, index, regimes(si, x[i], sj, x[j], mytree->time[node], mytree) ); 
							index++;
						}
					}
				}
			}
		}
	}
	free(x);
	free(node_begins);
}







/** Defines the transition matrices object used for the matrix-based recursion method */
void matrix_regimes(tree * mytree)
{
	/* Calc size of the data structure to hold all the matrices */
	size_t gridsize = mytree->gridsize;
	size_t n = mytree->n;
	size_t n_tips = (n+1)/2;
	size_t n_internal = (n-1)/2 - 1; // internals not including root
	size_t n_states = 2;
	size_t mem_block = n_states * (1+n_tips)*gridsize + n_internal* gsl_pow_2(gridsize)*gsl_pow_2(n_states);
	gsl_vector *regimes_matrices = gsl_vector_alloc(mem_block);

	printf("mem %d\n", mem_block);

	regimes_matrices_init(regimes_matrices, mytree);

	int * indices = indices_alloc_init(mytree);


//	gsl_vector_view V = gsl_vector_subvector(regimes_matrices, indices[2], gridsize*gridsize);
//	gsl_matrix_view M = gsl_matrix_view_array( (&V.vector)->data, gridsize, gridsize );
//	gsl_matrix_fprintf(stdout, &M.matrix, "%g");

	free(indices);  
	gsl_vector_free(regimes_matrices);

}

