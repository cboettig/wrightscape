/*  @file mvn.c
 *  @title multivariate normal random variables
 *  @author Carl Boettiger, <cboettig@gmail.com>
 *
 *  Based on the R function rmvnorm, from the mvtnorm package
 *  by Friedrich Leisch and Fabian Scheipl, implemented
 *  using the GSL libraries
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

int mvn(gsl_rng * rng, const gsl_vector * mean, gsl_matrix * covar, gsl_vector * ANS)
{
	int i;
	size_t n = mean->size;

	/* Calculate eigenvalues and eigenvectors of covar matrix */
	gsl_vector *eval = gsl_vector_alloc (n);
	gsl_matrix *evec = gsl_matrix_alloc (n, n);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
	gsl_eigen_symmv (covar, eval, evec, w);
	gsl_eigen_symmv_free (w);
//	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);


	
	/* Setup for: evec * matrix(diag(eval)) * transpose(evec)  */
	gsl_matrix *eval_mx = gsl_matrix_calloc (n, n);
	gsl_matrix * x_M = gsl_matrix_alloc (n,n);
	gsl_matrix * x_M_x = gsl_matrix_alloc (n,n);


	gsl_vector_view diagonal = gsl_matrix_diagonal(eval_mx);
	gsl_vector_memcpy(&diagonal.vector, eval);
	for(i=0;i<n;i++)
	{
		gsl_vector_set( &diagonal.vector, 
						i,  
						sqrt( gsl_vector_get(&diagonal.vector, i) )
					  );
	}



	/* evec * matrix(diag(eval)) * transpose(evec)  */
//	gsl_blas_dsymm (CblasLeft, CblasUpper, 
//					1.0, evec, eval_mx, 0.0, x_M);

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 
					1.0, evec, eval_mx, 0.0, x_M);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 
					1.0, x_M, evec, 0.0, x_M_x);

	
	gsl_matrix_free(x_M);
	gsl_matrix_free(eval_mx);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);

	gsl_vector * rnorms = gsl_vector_alloc(n);
	for(i=0;i<n;i++)
	{ 
		gsl_vector_set 
			( rnorms, i, 
			  gsl_ran_gaussian_ziggurat(rng, 1)
			);
	}

	gsl_blas_dgemv( CblasTrans, 1.0, x_M_x, rnorms, 0, ANS);
	gsl_vector_add(ANS, mean);
	gsl_matrix_free(x_M_x);
  gsl_vector_free(rnorms);

	return 0;
	/* answer provided through pass by reference */
}


/*
int main()
{
	const int n_tips = 2;
	gsl_rng * rng = gsl_rng_alloc(gsl_rng_default);
	gsl_vector * simdata = gsl_vector_alloc(n_tips);
	gsl_matrix * I = gsl_matrix_alloc(n_tips, n_tips);
	gsl_matrix_set_identity(I);
	gsl_vector * A = gsl_vector_calloc(n_tips);
	gsl_vector_set_all(A, 5.);
	mvn(rng, A, I, simdata);
	gsl_vector_fprintf(stdout, A, "%g");
	return 0;
}
*/

/* 
#define N 3
int main(void)
{

	double mean[] = {1,2,5};
	double sigma[] = {1, 0, 0,
					  0, 1, .9,
					  0, .9, 1};
	double out[N] = {0,0,0};

	gsl_vector_view mean_view = gsl_vector_view_array(mean, N);
	gsl_matrix_view sigma_view = gsl_matrix_view_array(sigma, N,N);
	gsl_vector_view out_view = gsl_vector_view_array(out, N);

	gsl_rng * rng = gsl_rng_alloc(gsl_rng_default);

	int i;
	for(i=0;i<10000;i++) {
		mvn(rng, &mean_view.vector, &sigma_view.matrix, &out_view.vector);
		printf("%g %g %g\n", out[0], out[1], out[2]);
	}
	return 0;
}
*/
