#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ilcplex/cplex.h>
#ifndef _MINIMIZATION_H_
#define _MINIMIZATION_H_

void minimization(contact *current);

static int quad_programming(double *x, double *H, double* P, double *C, double *B);

static int set_quad_problem(double *H, double* P, double *C, double *B, char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p,char **sense_p, int **matbeg_p, int **matcnt_p,int **matind_p, double **matval_p, double **lb_p,double **ub_p, int **qmatbeg_p, int **qmatcnt_p, int **qmatind_p, double **qmatval_p);

#endif
