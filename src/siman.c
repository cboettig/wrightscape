#include "tree.h"
#include <gsl/gsl_siman.h>

/* set up parameters for this simulated annealing run */
#define N_TRIES 200     	/* how many points do we try before stepping */
#define ITERS_FIXED_T 200	/* how many iterations for each T? */
#define STEP_SIZE 1. 		/* max step size in random walk */
#define K 1.0				/* Boltzmann constant */
#define T_INITIAL 0.01		/* initial temperature */
#define MU_T 1.003		    /* damping factor for temperature */
#define T_MIN .0078



  typedef struct {
	size_t size;
	double * data;
	void * s;
    } block;

  block * block_alloc(size_t n)
     {
		block * t = (block *) malloc(sizeof(block));
		t->data = (double *) malloc ((n+1) * sizeof(double));
		t->size = n;
		return t;
     }

  void block_free(block * t)
     {
		free(t->data);
		free(t);
     }

	void block_copy(void *inp, void *outp)
	{
		int i;
		block * in = (block *) inp;
		block * out = (block *) outp;

		for(i=0; i< in->size; i++){
			out->data[i] = in->data[i];
		}
		out->size = in->size;

		out->s = in-> s;
	}

	void * block_copy_construct(void *xp)
	{
		block * x = (block *) xp;
		block * y = block_alloc(x->size);
		block_copy(x, y);
		return y;
	}

	void block_destroy(void *xp)
	{
		block_free( (block *) xp);
	}


	/* distance function */
     double M1(void *xp, void *yp)
     {
       block * x = ((block *) xp);
       block * y = ((block *) yp);
	   int i;
	   double distance = 0;
	   for(i=0; i< x->size; i++){
        	distance += gsl_pow_2(x->data[i] - y->data[i]);
	   }
	   return sqrt(distance);
     }


	 double round(double);

	
	 /* Step function */ 
     void S1(const gsl_rng * r, void *xp, double step_size)
     {
       block *x = ((block *) xp);

	   int i = (int) round(gsl_rng_uniform(r)*x->size);

       x->data[i] = GSL_MAX(0, x->data[i]+gsl_ran_gaussian_ziggurat(r,step_size)); //gsl_rng_uniform(r) * 2 * step_size - step_size;
//	   printf("%.2g\n", x->data[i]);

     }

	 /* print function */
     void P1(void *xp)
     {
		block * x = (block *) xp;
		int i;
		for(i=0;i< x->size; i++){
			printf("%6.2lf ", x->data[i]);
		}
     }


	/* distance function */
     double dist(void *xp, void *yp)
     {
       tree * x = ((tree *) xp);
       tree * y = ((tree *) yp);
	   int i;
	   double distance = 0;
	   for(i=0; i< x->npars; i++){
			if(x->fitpar[i]==1){
				distance += gsl_pow_2(x->pars[i] - y->pars[i]);
			}
	   }
	   return sqrt(distance);
     }


	 double round(double);


	double E1(void * x){
		block * myblock = (block *) x; 
		tree * mytree = (tree *) (myblock->s);

		/* Write block parameters into mytree  */
		int i, j=0;
		for(i=0; i<mytree->npars; i++){
			if(mytree->fitpar[i] == 1){
				mytree->pars[i] = myblock->data[j];
				if(mytree->pars[i] < 0) return 1000; // Force parameters to be positive only!
				j++;
			}
		}

		return -ou_likelihood(mytree);
	}



	void siman(tree * mytree, gsl_rng * rng) 
     {
		block * myblock = block_alloc(mytree->nfreepars);

		/* initialize block */
		int i, j = 0;
		for(i = 0; i < mytree->npars; i++){
			if(mytree->fitpar[i] == 1){ 
				myblock->data[j] = mytree->pars[i];
				j++;
			}
		}
		myblock->s = mytree;

		gsl_siman_params_t siman_params
		   = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
			  K, T_INITIAL, MU_T, T_MIN};

       gsl_siman_solve(rng, myblock, E1, S1, M1, P1, block_copy, block_copy_construct, block_destroy,0, siman_params);
   	   P1(myblock); printf("\n");

     }

