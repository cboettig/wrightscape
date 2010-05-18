#include "tree.h"		

/* The vector for the search space,then any other params*/
double my_f (const gsl_vector *v, void *params)
{
	tree * mytree = (tree *) params;
	int i, j=0;
	for(i = 0; i < mytree->npars; i++){
		if(mytree->fitpar[i] == 1){ 
			mytree->pars[i] = gsl_vector_get(v,j);
			if(mytree->pars[i] < 0) return 1000; // Force parameters to be positive only!
			j++;
		}
	}
	return -ou_likelihood(mytree);
}



void multimin(void * params)
{
	tree * mytree = (tree *) params;

	size_t iter = 0;
	int status;
	double size;
	int n_pars= 0;
	int i; 
	for(i=0; i<mytree->npars; i++){
		n_pars += mytree->fitpar[i];
	}
	
	/* Declare minimizer type and allocate it */
	gsl_multimin_fminimizer * s = 
		gsl_multimin_fminimizer_alloc
		(
			gsl_multimin_fminimizer_nmsimplex2, n_pars
		);

	/* Set initial step sizes to a vector of INIT_STEP size) */
	gsl_vector *ss = gsl_vector_alloc (n_pars);
	gsl_vector_set_all (ss, INIT_STEP);


	/* Extract the parameters vector and write into a gsl_vec to search over */

	gsl_vector *x = gsl_vector_alloc (n_pars);
	for(i=0; i < n_pars; i++){
		if(mytree->fitpar[i] == 1){
			gsl_vector_set(
				x,
				i, 
				mytree->pars[i]);
		}
	}


	/* Initialize method and iterate */
	gsl_multimin_function minimize_this;
	minimize_this.n = n_pars;  // dimension
	minimize_this.f = &my_f;
	minimize_this.params = params;

	gsl_multimin_fminimizer_set (s, &minimize_this, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, ERR_TOL);
		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}
			
	} while (status == GSL_CONTINUE && iter < MAX_ITER);

	/* Print out status, whether or not it converged */
	size = gsl_multimin_fminimizer_size (s);
	printf ("iter: %5d, -llik = %3.3f, size = %.3f, par values = ",
		   iter,
		   s->fval, 
		   size );
	for(i = 0; i < n_pars; i++)
		printf (" %1lf, ", gsl_vector_get (s->x, i) );
	printf("\n");

   /* Write the values of the gsl_vec x back into the parameter list */
	for(i=0; i < n_pars; i++){
		if(mytree->fitpar[i] == 1){
			gsl_vector_set(
				x,
				i, 
				mytree->pars[i]);
		}
	}

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

}

