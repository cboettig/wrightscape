/* Limits of integration*/

/* NB: A wide grid range is much more important than a fine grid! */
#define LLIM -12.0
#define ULIM 15.0
#define GRIDSIZE 51

/* multimin parameters */
#define INIT_STEP .1
#define MAX_ITER 1500
#define ERR_TOL 1e-6


#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


struct treestruct {
		int n;
		double * time;
		double * trait;
		int * state;	// regime indicator for each node
		int * ancestor; // sufficent to define the tree topology
		// For convenience/speed, can be computed from ancestor vector 
		int * left;	
		int * right;	
		int * tip;		

		size_t nstates;	// number of regimes
		
		size_t npars;		
		double * pars;
		int * fitpar;
		size_t nfreepars;

		size_t gridsize;
};

typedef struct treestruct tree;

/* functions in tree.c */
tree * tree_alloc(const size_t n, const size_t npars);
void tree_init(tree *t, const double * times, const int * ancestors, 
	const double * traits, const int * states, const int nstates, const int n, 
	const double *pars, const int * fitpars, const int gridsize);
void * tree_copy_construct(tree *in);
void tree_copy(tree *in, tree *out);
void tree_print(tree *t);
void tree_free(tree * t);

/* functions in multimin.c*/
double my_f (const gsl_vector *v, void *params);
double multimin(void * params);

/* functions in models.c*/
void matrix_regimes(tree * mytree);

/* (external) functions in siman */
double siman(tree * mytree, gsl_rng * rng);


/** Defines the transition matrices object used for the matrix-based recursion method */
typedef struct {
	size_t n;
	size_t gridsize;
	double ** matrix; // list of arrays
	double ** sofar;

	// we'll reuse these product vectors for each node
	double * product;
	double * grid;
} matrices;

/* functions in matrix_method.c */
double matrix_likelihood(tree * mytree);

/* functions in linearsoln.c */
double bm_likelihood(tree * mytree);
double ou_likelihood(tree * mytree);
double bm_ancestral_likelihood(tree * mytree);

