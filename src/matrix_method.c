#include "tree.h"

/** Computes transition density directly from inputs without reference to the tree
 * Eventually this should really decide which model is appropriate, not hardwire it.   
 *
 * */


//double w_direct(double x, double y, double t, double * params){ return exp(-gsl_pow_2(x-y)/(2*t*params[0]) )/sqrt(2*M_PI*t*params[0]); }

 
double w_direct(double x, double y, double t, double * params)
{ 
	double sigma2 = params[0];
	double alpha = params[1];
	double theta = params[2];
	double mean = x*exp(-alpha*t) + theta*(1-exp(-alpha*t)) ;
	double var = sigma2/(2*alpha)*(1-exp(-2*alpha*t) );
	return exp(-gsl_pow_2(y-mean)/(2*var) )/sqrt(2*M_PI*var); 
}


matrices * matrices_alloc(tree * mytree)
{
	size_t n = mytree->n, gridsize = mytree->gridsize;
	matrices * mymatrices = (matrices *) malloc(sizeof(matrices));
	mymatrices->matrix = (double **) malloc(gridsize * sizeof(double *));  
	mymatrices->sofar = (double **) malloc(gridsize * sizeof(double *));  
	int i;

	// first do the root
	mymatrices->matrix[0] = (double *) malloc(gridsize * sizeof(double)); // root gets column vector
	mymatrices->sofar[0] = (double *) malloc(gridsize * sizeof(double));
	for(i=1;i<n;i++){ //start after root
		mymatrices->sofar[i] = (double *) malloc(gridsize * sizeof(double));
		if(mytree->tip[i]){ // tips get column vectors
			mymatrices->matrix[i] = (double *) malloc(gridsize * sizeof(double));
		}
		else { // internals get matrices
			mymatrices->matrix[i] = (double *) malloc((gridsize*gridsize) * sizeof(double));
		}
	}

	mymatrices->product = (double *) malloc(gridsize * sizeof(double));  // products are row vectors
	mymatrices->grid = (double *) malloc(gridsize * sizeof(double));
	mymatrices->gridsize = gridsize;
	mymatrices->n = n;
	return mymatrices;
}

void matrices_free(matrices * mymatrices)
{
	int i;
	for(i=0;i<mymatrices->n;i++){
		free(mymatrices->matrix[i]);
		free(mymatrices->sofar[i]);
	}
	free(mymatrices->matrix);
	free(mymatrices->sofar);
	free(mymatrices->product);
	free(mymatrices->grid);
	free(mymatrices);
}

void matrices_init(tree * mytree, matrices * mymatrices)
{
	int i,j,k, s=0;
	double x = mytree->trait[0];
	size_t gridsize = mymatrices->gridsize;
	double delta_grid = (double) (ULIM-LLIM)/(gridsize-1);
	for(i=0;i<gridsize;i++){
		mymatrices->product[i] = 0;
		mymatrices->grid[i] = (double) delta_grid*i + LLIM;
	}

	// root from value in tree
	s = 1;
	for(i=0;i<mymatrices->gridsize;i++){
		(mymatrices->matrix[0])[i] = 0;
		if(i*delta_grid+LLIM > x && s){
			(mymatrices->matrix[0])[i] = 1;
			s = 0;
		}
	}


	for(k=1;k<mytree->n;k++){ //skip root
		if(mytree->tip[k]){ // tips get column vectors
			for(i=0;i<gridsize;i++){
				(mymatrices->matrix[k])[i] = 
					w_direct (	mymatrices->grid[i], 
								mytree->trait[k], 
								mytree->time[k], 
								mytree->pars) ; 
			}
		} else {
			for(i=0;i<gridsize;i++){
				for(j=0;j<gridsize;j++){
					(mymatrices->matrix[k])[i*gridsize+j] = 
						w_direct (	mymatrices->grid[i],
									mymatrices->grid[j],
									mytree->time[k], 
									mytree->pars) * delta_grid; 
				}
			}
		}
	}
}


// Consider generating matrices by diagonalizing the elemental process and then raising to power t
// This would enable sparse matrix multiply
//
// Call on root using: matrix_recurse(mymatrices->matrix[0], 0, mytree, mymatrices), this way gets root sofar as row vector

/*
 *	
 *	g
	|\
	f \
	|\ \
	e \ \
	|\ \ \
    a b c d

*/

double matrix_recurse(int i, tree * mytree, matrices * mymatrices)
{
	int j;
	int left, right;
	double ans;
	size_t gridsize = mymatrices->gridsize;

	if(mytree->tip[i]){
		for(j=0;j<gridsize;j++){
			(mymatrices->sofar[i])[j] = (mymatrices->matrix[i])[j]; // this could be done by init, and avoid alloc matrix[i] for tips 
		}
		return 0;
	} else {
		/* Recurse on both the left and the right */
		left = mytree->left[i];
		right = mytree->right[i];
		matrix_recurse(left, mytree, mymatrices); 
		matrix_recurse(right, mytree, mymatrices);


		/* elementwise multiply the left and right, which are tips or reduced by integral to being a tip */
		gsl_vector_view LEFT = gsl_vector_view_array(mymatrices->sofar[left], gridsize);
		gsl_vector_view RIGHT = gsl_vector_view_array(mymatrices->sofar[right], gridsize);
		gsl_vector_mul(&LEFT.vector, &RIGHT.vector);

		if(i==0){
			/* If we're finally back to the root, return the dot product of the known root with the vector */
			gsl_vector_view R = gsl_vector_view_array( mymatrices->matrix[0], gridsize) ;
			gsl_blas_ddot(&R.vector, &LEFT.vector, &ans);
			return log(ans); 

		} else {
			/* Integral is the vector product of the elementwise product of tips with matrix of the ancestor */
			gsl_matrix_view M = gsl_matrix_view_array( mymatrices->matrix[i], gridsize, gridsize ) ;
			gsl_vector_view ANS = gsl_vector_view_array(mymatrices->sofar[i], gridsize);
			gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &LEFT.vector, 0, &ANS.vector);
			/* Copy the product into the first row of the matrix for that set, avoids using seperate vector */

		}
	}
	return 1;
}



/* Return log likelihood using matrix method.  Commented out code allows integrating over root */
double matrix_likelihood(tree * mytree)
{
	matrices *mymatrices = matrices_alloc(mytree);
	matrices_init(mytree, mymatrices);

	/* likelihood given a fixed root parameter */
	double ans = matrix_recurse(0, mytree, mymatrices) ;

	matrices_free(mymatrices);	
	return ans;

	/* likelihood integrating over each possible value for root */
/*	
 	double prob_over_all_roots = 0;
	int i,j;

	for(i=0;i<mymatrices->gridsize;i++){
		for(j=0;j<mymatrices->gridsize;j++){
				(mymatrices->sofar[0])[j] = 0;
		}
		(mymatrices->sofar[0])[i] = 1;
			prob_over_all_roots += exp(matrix_recurse(0, mytree, mymatrices)-log(gridsize)) ;  //  e^r/gridsize
	}
	matrices_free(mymatrices);	
	return log(prob_over_all_roots);
*/

}







void matrices_print(tree * mytree, matrices * mymatrices)
{
	int i,j,k, s=0;
	size_t gridsize = mymatrices->gridsize;
	printf("grid:\n");
	for(i=0;i<gridsize;i++){
		printf("%lf ", mymatrices->grid[i]);
	}
	printf("\n\n");

	for(k=1;k<mytree->n;k++){ //skip root
		s=0;
		if(mytree->tip[k]){ // tips get column vectors
			for(i=0;i<gridsize;i++){
				printf("%lf, ", (mymatrices->matrix[k])[s] );
				s++;
			}
			printf("\n");
		} else {
			for(i=0;i<gridsize;i++){
				for(j=0;j<gridsize;j++){
					printf("%lf, ", (mymatrices->matrix[k])[s] );
					s++;
				}
				printf("\n");
			}
			printf("\n");
		}
		printf("\n");
	}
}


