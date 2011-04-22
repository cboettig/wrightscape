/**
* @file likelihood.c
* @brief base functions to calculate the likelihood for general multitype OU processes on a phylogenetic tree
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include "optimizers.h"
#include "mvn.h"


typedef struct {
	size_t n_nodes;
	size_t n_regimes;
	int * ancestor;
	int * regimes;
	double * branch_length;
	double * traits;
	double * alpha;
	double * theta;
	double * sigma;
	double * Xo;
	int * lca_matrix;
} tree;


int get_lca (int i, int j, int n_nodes, const int * ancestor, 
             const double * branch_length, double * sep);
void simulate (const gsl_rng * rng, tree * mytree);
double calc_lik (const double *Xo, const double * alpha, const double * theta,
                 const double * sigma, const int * regimes, 
                 const int * ancestor, const double * branch_length, 
                 const double * traits, int n_nodes, int * lca_matrix);


