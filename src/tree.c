#include "tree.h"
/* functions for handling the tree data */
tree * tree_alloc(const size_t n, const size_t npars)
{
		tree * t = (tree   *) malloc(sizeof(tree));
		t->time  = (double *) malloc ((n+1) * sizeof(double));
		t->trait = (double *) malloc ((n+1) * sizeof(double));
		t->state = (int *) malloc ((n+1) * sizeof(int));
		t->ancestor = (int *) malloc ((n+1) * sizeof(int));
		t->left = (int *) malloc ((n+1) * sizeof(int));
		t->right = (int *) malloc ((n+1) * sizeof(int));
		t->tip = (int *) malloc ((n+1) * sizeof(int));
		t->n = n;

		t->pars = (double *) malloc ((npars+1) * sizeof(double));
		t->fitpar = (int *) malloc ((npars+1) * sizeof(int));
		t->npars = npars;

		return t;
}

void tree_free(tree * t)
{
	free(t->time);
	free(t->trait);
	free(t->state);
	free(t->ancestor);
	free(t->left);
	free(t->right);
	free(t->tip);
	free(t->pars);
	free(t->fitpar);
	free(t);
}

/* Ancestor information kept by binary indexing  */ 
void tree_init(tree *t, const double * times, const int * ancestors, const double * traits, const int * states, const int nstates, const int n, const double * pars, const int * fitpars, const int gridsize)
{
	int i, j=0, left_empty=1;
	for(i=0; i < n; i++){ 
		left_empty = 1;
		t->time[i] = times[i];
		t->trait[i] = traits[i];
		t->state[i] = states[i]; 
		t->ancestor[i] = ancestors[i];
		t->left[i] = 0;
		t->right[i] = 0;
		t->tip[i] = 0;
		for(j=0; j < n; j++){
			if(i == ancestors[j] ){
				if(left_empty){ 
					t->left[i] = j;
					left_empty = 0;
				} else { 
					t->right[i] = j;
				}
			}
		}
		if (t->left[i] == 0) t->tip[i] = 1;
	}
	t->nstates = nstates;
	t->nfreepars = 0;
	j=0;
	for(i=0; i < t->npars; i++){
		t->pars[i] = pars[i];
		t->fitpar[i] = fitpars[i];
		t->nfreepars += fitpars[i];
	}
	t->gridsize = gridsize;
}


void tree_print(tree * t)
{
	int i;
	for(i=0; i< t->n; i++){
			printf("%d %.2lf %.2lf \t %4d %4d %4d  %4d\n", i, t->time[i], t->trait[i], t->ancestor[i], t->left[i], t->right[i], t->tip[i] );
	}
}





void tree_copy(tree * in, tree * out)
{
	int i;	
	for(i=0; i<= in->n ; i++){
		out->time[i]  =  in->time[i] ;
		out->trait[i] = in->trait[i] ;
		out->state[i] = in->state[i] ;
		out->ancestor[i] = in->ancestor[i] ;
		out->left[i] = in->left[i] ;
		out->right[i] = in->right[i] ;
		out->tip[i] = in->tip[i] ;
	}
	for(i=0; i<= in->npars ; i++){
		out->pars[i]  =  in->pars[i] ;
		out->fitpar[i]  =  in->fitpar[i] ;
	}	
	out->nstates = in->nstates;
}

void * tree_copy_construct(tree * in)
{
	tree * t = tree_alloc(in->n, in->npars);
	int i;	
	for(i=0; i<= in->n ; i++){
		t->time[i]  =  in->time[i] ;
		t->trait[i] = in->trait[i] ;
		t->state[i] = in->state[i] ;
		t->ancestor[i] = in->ancestor[i] ;
		t->left[i] = in->left[i] ;
		t->right[i] = in->right[i] ;
		t->tip[i] = in->tip[i] ;
	}
	for(i=0; i<= in->npars ; i++){
		t->pars[i]  =  in->pars[i] ;
		t->fitpar[i]  =  in->fitpar[i] ;
	}	
	t->nstates = in->nstates;
	return t;
}


