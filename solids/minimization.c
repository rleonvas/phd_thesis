#include "functions.h"
#include "minimization.h"

void minimization(contact *current){
/*
	Minimize F(x) = 0.5*x'(2M + PHI)x - x'(2MU + Pext) = x'Hx - x'P
	subject to
		C.U'/2 - C.x' <= 0  :  C.x' >= C.U'/2 : C.x' >= B
*/

	int i, j, k, l, ind_one, ind_two, status;
	double *x, *H = NULL, *C = NULL, *P = NULL, *B = NULL, *val_aux = NULL;
	double value, r[3][12], s[3][12], t[3][12];	
	contact *aux_a,*aux_b,*aux_c;

	/*printf("\n");
	for(i = 0; i < N; i++){
		printf("SOLID %d: ",i);
		for(j = 0; j < solid[i].cantBlock; j++)
			printf("%d ",solid[i].idBlock[j]);
		printf(" -- CB : %d\n",solid[i].cantBlock);
	}*/

	val_aux = (double*)calloc(36,sizeof(double));
	//Computing mass and inertia matrix M as a vector.
	H = (double*)calloc(36*N*N,sizeof(double));
	for(i = 0;i < N;i++){
		for(j = 0;j < 3;j++){
			*(H + (6*i + j)*(6*N + 1)) = 2.*solid[i].mass;
			for(k = 0;k < 3;k++)
				*(H + (6*i + j + 3)*(6*N + 1) - j + k) = 2.*solid[i].inertia_matrix[j][k];
		}
	}

	/*printf("MATRIZ M\n");
	displayMatrix(M, 6*N, 6*N);*/

	//Computing vector --> P = 2MU + Pext (linear part)
	P = (double*)calloc(6*N,sizeof(double));
	for(i = 0;i < N;i++){
		for(j = 0;j < 3;j++){
			*(P + 6*i + j) += 2.*solid[i].mass*solid[i].vel[j];
			value = 0.;
			for(k = 0;k < 3;k++){
				value += solid[i].inertia_matrix[j][k]*solid[i].omg[k];
			}
			*(P + 6*i + j + 3) += 2.*value;
		//Adding Pext
		*(P + 6*i + 2) -= solid[i].mass*G*timestep;
		}
	}
	/*printf("VECTOR P\n");
	displayMatrix(P, 6*N, 1);*/

  aux_a = current;
	aux_b = current;
	aux_c = current;
	//Computing the constraint matrix --> vector C
	C = (double*)calloc(6*nc*N,sizeof(double));
	for(i = 0; i < nc;i++){
		if (current->sol_two < 0){
			for(k = 0; k < 3; k++)
				*(C + 6*N*i + 6*(current->sol_one) + k) = current->normal_vec[k];
*(C + 6*N*i + 6*(current->sol_one) + 3) = current->normal_vec[2]*current->gp_one[1] - current->normal_vec[1]*current->gp_one[2];
*(C + 6*N*i + 6*(current->sol_one) + 4) = current->normal_vec[0]*current->gp_one[2] - current->normal_vec[2]*current->gp_one[0];
*(C + 6*N*i + 6*(current->sol_one) + 5) = current->normal_vec[1]*current->gp_one[0] - current->normal_vec[0]*current->gp_one[1];
		}
		else{
			for(k = 0; k < 3; k++){
				*(C + 6*N*i + 6*(current->sol_one) + k) = current->normal_vec[k];
				*(C + 6*N*i + 6*(current->sol_two) + k) = -(current->normal_vec[k]);
			}
*(C + 6*N*i + 6*(current->sol_one) + 3) = current->normal_vec[2]*current->gp_one[1] - current->normal_vec[1]*current->gp_one[2];
*(C + 6*N*i + 6*(current->sol_one) + 4) = current->normal_vec[0]*current->gp_one[2] - current->normal_vec[2]*current->gp_one[0];
*(C + 6*N*i + 6*(current->sol_one) + 5) = current->normal_vec[1]*current->gp_one[0] - current->normal_vec[0]*current->gp_one[1];
*(C + 6*N*i + 6*(current->sol_two) + 3) = current->normal_vec[1]*current->gp_two[2] - current->normal_vec[2]*current->gp_two[1];
*(C + 6*N*i + 6*(current->sol_two) + 4) = current->normal_vec[2]*current->gp_two[0] - current->normal_vec[0]*current->gp_two[2];
*(C + 6*N*i + 6*(current->sol_two) + 5) = current->normal_vec[0]*current->gp_two[1] - current->normal_vec[1]*current->gp_two[0];
		}
		current = current->next;
	}
	/*printf("MATRIZ C\n");
	displayMatrix(C, nc, 6*N);*/

	/*Computing right side of constraints --> vector B = C.U'/2 */
	B = (double*)calloc(nc,sizeof(double));
	for(i = 0;i < nc;i++){
		if (aux_a->sol_two < 0){
			for(k = 0; k < 3; k++){
				*(B + i) += C[6*N*i + 6*(aux_a->sol_one) + k]*solid[aux_a->sol_one].vel[k];
				*(B + i) += C[6*N*i + 6*(aux_a->sol_one) + k + 3]*solid[aux_a->sol_one].omg[k];
			}
		}
		else{
			for(k = 0; k < 3; k++){
				*(B + i) += C[6*N*i + 6*(aux_a->sol_one) + k]*solid[aux_a->sol_one].vel[k] + C[6*N*i + 6*(aux_a->sol_two) + k]*solid[aux_a->sol_two].vel[k];
				*(B + i) += C[6*N*i + 6*(aux_a->sol_one) + k + 3]*solid[aux_a->sol_one].omg[k] + C[6*N*i + 6*(aux_a->sol_two) + k + 3]*solid[aux_a->sol_two].omg[k];
			}
		}
    *(B + i) = .5*B[i];
    aux_a = aux_a->next;
	}
	/*printf("VECTOR B\n");
	displayMatrix(B, nc, 1);*/

	//------------------------------------ H += kn*PHI_N ------------------------------------
	for(i = 0;i < nc;i++){
		if(aux_b->sol_two < 0){
			for(j = 0;j < 6; j++){
				for(k = 0;k < 6; k++){
*(H + 6*(aux_b->sol_one)*(6*N + 1) + 6*N*j + k) += kn*C[6*N*i + 6*(aux_b->sol_one) + j]*C[6*N*i + 6*(aux_b->sol_one) + k];
				}
			}
		}
		else{
      for(j = 0;j < 6; j++){
				for(k = 0;k < 6; k++){
					*(H + 6*(aux_b->sol_one)*(6*N + 1) + 6*N*j + k) += kn*C[6*N*i + 6*(aux_b->sol_one) + j]*C[6*N*i + 6*(aux_b->sol_one) + k];
					*(H + 6*(aux_b->sol_two)*(6*N + 1) + 6*N*j + k) += kn*C[6*N*i + 6*(aux_b->sol_two) + j]*C[6*N*i + 6*(aux_b->sol_two) + k];
					*(val_aux + 6*j + k) = kn*C[6*N*i + 6*(aux_b->sol_one) + j]*C[6*N*i + 6*(aux_b->sol_two) + k];
				}
			}
			for(j = 0;j < 6; j++){
				for(k = 0;k < 6; k++){
					*(H + 36*N*(aux_b->sol_one) + 6*(aux_b->sol_two) + 6*N*j + k) += val_aux[6*j + k];
					*(H + 36*N*(aux_b->sol_two) + 6*(aux_b->sol_one) + 6*N*j + k) += val_aux[6*k + j];
         }
			}
		}
	  aux_b = aux_b->next;
	}

	//------------------------------------ H += kt*PHI_T ------------------------------------
	for(i = 0;i < nc;i++){
		if(aux_c->sol_two < 0){
      /*Computing matrices R and S.*/
      for(j = 0;j < 3;j++){
		    for(k = 0;k < 6;k++){
          r[j][k] = (j == k) ? 1. : 0.;
          s[j][k] = C[6*N*i + 6*(aux_c->sol_one) + k]*aux_c->normal_vec[j];
		    }
      }
      r[0][4] = aux_c->gp_one[2];
      r[0][5] = -(aux_c->gp_one[1]);
      r[1][3] = -(aux_c->gp_one[2]);
      r[1][5] = aux_c->gp_one[0];
      r[2][3] = aux_c->gp_one[1];
      r[2][4] = -(aux_c->gp_one[0]);

			/*Computing matrix T = R - S.*/
			for(j = 0;j < 3;j++){
       for(k = 0;k < 6;k++)
				t[j][k] = r[j][k] - s[j][k];
      }
			/*Computing T'*T*/
			for(j = 0;j < 6;j++){
				for(k = 0;k < 6;k++){
          for(l = 0;l < 3;l++)
            *(H + 6*(aux_c->sol_one)*(6*N + 1) + 6*N*j + k) += kt*t[l][j]*t[l][k];
  				}
			}
		}
    else{
      /*Computing matrices R and S.*/
      for(j = 0;j < 3;j++){
        for(k = 0;k < 6;k++){
          r[j][k] = (j == k) ? 1. : 0.;
          r[j][k + 6] = (j == k) ? -1. : 0.;
          s[j][k] = C[6*N*i + 6*(aux_c->sol_one) + k]*aux_c->normal_vec[j];
          s[j][k + 6] = C[6*N*i + 6*(aux_c->sol_two) + k]*aux_c->normal_vec[j];
        }
      }
      r[0][4] = aux_c->gp_one[2];
      r[0][5] = -(aux_c->gp_one[1]);
      r[1][3] = -(aux_c->gp_one[2]);
      r[1][5] = aux_c->gp_one[0];
      r[2][3] = aux_c->gp_one[1];
      r[2][4] = -(aux_c->gp_one[0]);
      r[0][10] = -(aux_c->gp_two[2]);
      r[0][11] = aux_c->gp_two[1];
      r[1][9] = aux_c->gp_two[2];
      r[1][11] = -(aux_c->gp_two[0]);
      r[2][9] = -(aux_c->gp_two[1]);
      r[2][10] = aux_c->gp_two[0];

			/*Computing matrix T = R - S.*/
			for(j = 0;j < 3;j++){
       for(k = 0;k < 12;k++)
		    t[j][k] = r[j][k] - s[j][k];
      }

			/*Computing T'*T*/
			for(j = 0;j < 6;j++){
				for(k = 0;k < 6;k++){
					val_aux[6*j + k] = 0.;
          for(l = 0;l < 3;l++){
		      	*(H + 6*(aux_c->sol_one)*(6*N + 1) + 6*N*j + k) += kt*t[l][j]*t[l][k];
		        *(H + 6*(aux_c->sol_two)*(6*N + 1) + 6*N*j + k) += kt*t[l][j + 6]*t[l][k + 6];
						*(val_aux + 6*j + k) += kt*t[l][j]*t[l][k + 6];		        
          }
  			}
			}
			for(j = 0;j < 6; j++){
				for(k = 0;k < 6; k++){
					*(H + 36*N*(aux_c->sol_one) + 6*(aux_c->sol_two) + 6*N*j + k) += val_aux[6*j + k];
					*(H + 36*N*(aux_c->sol_two) + 6*(aux_c->sol_one) + 6*N*j + k) += val_aux[6*k + j];
         }
			}
   	}
	  aux_c = aux_c->next;
	}	

	/*H must be symmetric*/
	for(i = 0; i < N; i++){
		for(j = 0; j < solid[i].cantBlock; j++){
			ind_one = 36*N*i + 6*solid[i].idBlock[j];
			ind_two = 36*N*solid[i].idBlock[j] + 6*i;
			for(l = 0; l < 6; l++){
				for(k = 0; k < 6; k++){
					H[ind_one + 6*N*l + k] = .5*(H[ind_one + 6*N*l + k] + H[ind_two + 6*N*k + l]);
		    	H[ind_two + 6*N*k + l] = H[ind_one + 6*N*l + k];
				}
			}
		}
	}

	/*printf("MATRIZ H\n");
	displayMatrix(H, 6*N, 6*N);*/

	x = (double*)malloc(6*N*sizeof(double));
	status = quad_programming(x, H, P, C, B);

	/*printf("VECTOR X\n");
	displayMatrix(x, 6*N, 1);*/

	//Updating velocities with contacts.
	for(i = 0;i < N;i++){
		for(j = 0;j < 3;j++){
			solid[i].vel[j] = 2.*x[6*i + j] - solid[i].vel[j];
			solid[i].omg[j] = 2.*x[6*i + j + 3] - solid[i].omg[j];
		}
	}

	//Freeing memory.
	free_and_null((char **) &val_aux);
	free_and_null((char **) &H);
	free_and_null((char **) &C);
	free_and_null((char **) &P);
	free_and_null((char **) &B);	
	free_and_null((char **) &current);
	free_and_null((char **) &aux_a);
	free_and_null((char **) &aux_b);
	free_and_null((char **) &aux_c);
	free_and_null((char **) &x);		
}

static int quad_programming(double *x, double *H, double* P, double *C, double *B){

	// Declare pointers for the variables and arrays that will contain the data which define the LP problem.
	// The setproblemdata() routine allocates space for the problem data.

	char     *probname = NULL;
	int      numcols;
	int      numrows;
	int      objsen;
	double   *obj = NULL;
	double   *rhs = NULL;
	char     *sense = NULL;
	int      *matbeg = NULL;
	int      *matcnt = NULL;
	int      *matind = NULL;
	double   *matval = NULL;
	double   *lb = NULL;
	double   *ub = NULL;
	int      *qmatbeg = NULL;
	int      *qmatcnt = NULL;
	int      *qmatind = NULL;
	double   *qmatval = NULL;

	// Declare and allocate space for the variables and arrays where we will store the optimization results including
	// the status, objective value, variable values, dual values, row slacks and variable reduced costs.

	int 		solstat;
	double		objval;
	double		pi[nc];
	double		slack[nc];
	double		dj[6*N];

	CPXENVptr     env = NULL;
	CPXLPptr      lp = NULL;
	int           status;
	int           i, j;
	int           cur_numrows, cur_numcols;

	// Initialize the CPLEX environment
	env = CPXopenCPLEX(&status);
	// If an error occurs, the status value indicates the reason for failure.
	// A call to CPXgeterrorstring will produce the text of the error message.
	// Note that CPXopenCPLEX produces no output, so the only way to see the cause of the error is to use
	// CPXgeterrorstring.  For other CPLEX routines, the errors will be seen if the CPX_PARAM_SCRIND indicator
	// is set to CPX_ON.

	if(env == NULL){
		char errmsg[1024];
		fprintf(stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring(env, status, errmsg);
		fprintf(stderr, "%s", errmsg);
		goto TERMINATE;
	}

	// Turn on output to the screen
	status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	if(status){
		fprintf(stderr, "Failure to turn on screen indicator, error %d.\n", status);
		goto TERMINATE;
	}

	// Fill in the data for the problem.
	status = set_quad_problem(H, P, C, B, &probname, &numcols, &numrows, &objsen, &obj, &rhs, &sense, &matbeg, &matcnt, &matind, &matval, &lb, &ub, &qmatbeg, &qmatcnt, &qmatind, &qmatval);

	if(status){
		fprintf(stderr, "Failed to build problem data arrays.\n");
		goto TERMINATE;
	}

	// Create the problem.
	lp = CPXcreateprob(env, &status, probname);
	// A returned pointer of NULL may mean that not enough memory was available or there was some other problem.
	// In the case of failure, an error message will have been written to the error channel from inside CPLEX.
	// The setting of the parameter CPX_PARAM_SCRIND causes the error message to appear on stdout.
	if(lp == NULL) {
		fprintf(stderr, "Failed to create problem.\n");
		goto TERMINATE;
	}

	// Now copy the LP part of the problem data into the lp
	status = CPXcopylp(env, lp, numcols, numrows, objsen, obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
	if(status){
		fprintf(stderr, "Failed to copy problem data.\n");
		goto TERMINATE;
	}

	// Now copy the QUAD part of the problem data into the lp
	status = CPXcopyquad(env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
	if(status){
		fprintf (stderr, "Failed to copy quadratic matrix.\n");
		goto TERMINATE;
	}

  status = CPXsetintparam(env, CPX_PARAM_QPMETHOD, 4);
	if(status){
		fprintf(stderr, "\nFailed to set the quadratic method.", status);
		goto TERMINATE;
	}

	status = CPXsetintparam(env, CPX_PARAM_BARALG, 3);
	if(status){
		fprintf(stderr, "\nFailed to set the barrier algorithm.", status);
		goto TERMINATE;
	}

  status = CPXsetdblparam(env, CPX_PARAM_BAREPCOMP, 1e-3);
	if(status){
		fprintf(stderr, "\nFailed to set the convergence parameter.", status);
		goto TERMINATE;
	}

	// Optimize the problem and obtain solution.
	status = CPXqpopt(env, lp);
	if(status){
		fprintf(stderr, "Failed to optimize QP.\n");
		goto TERMINATE;
	}

	status = CPXsolution(env, lp, &solstat, &objval, x, pi, slack, dj);
	if(status){
		fprintf(stderr, "Failed to obtain solution.\n");
		goto TERMINATE;
	}

	// Write the output to the screen.
	//printf("\nSolution status = %d\n", solstat);
 	//printf("Solution value  = %f\n\n", objval);

	// The size of the problem should be obtained by asking CPLEX what the actual size is,
	// rather than using what was passed to CPXcopylp.
	// cur_numrows and cur_numcols store the current number of rows and columns, respectively.

	//cur_numrows = CPXgetnumrows(env, lp);
	//cur_numcols = CPXgetnumcols(env, lp);
	//for(i = 0; i < cur_numrows; i++)
	//	printf("Row %d:  Slack = %.10f  Pi = %.10f\n", i, slack[i], pi[i]);

	//for(j = 0; j < cur_numcols; j++)
	//	printf("Column %d:  Value = %.10f  Reduced cost = %.10f\n",j, x[j], dj[j]);

	// Finally, write a copy of the problem to a file.
	//status = CPXwriteprob(env, lp, "quad_prog.lp", NULL);
	//if(status){
	//	fprintf(stderr, "Failed to write LP to disk.\n");
	//	goto TERMINATE;
	//}

TERMINATE:
        // Free up the problem as allocated by CPXcreateprob, if necessary
	if(lp != NULL) {
		status = CPXfreeprob(env, &lp);
		if(status){
	 		fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	// Free up the CPLEX environment, if necessary
	if (env != NULL){
		status = CPXcloseCPLEX(&env);
		if (status){
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}
	// Free up the problem data arrays, if necessary.
	free_and_null((char **) &probname);
	free_and_null((char **) &obj);
	free_and_null((char **) &rhs);
	free_and_null((char **) &sense);
	free_and_null((char **) &matbeg);
	free_and_null((char **) &matcnt);
	free_and_null((char **) &matind);
	free_and_null((char **) &matval);
	free_and_null((char **) &lb);
	free_and_null((char **) &ub);
	free_and_null((char **) &qmatbeg);
	free_and_null((char **) &qmatcnt);
	free_and_null((char **) &qmatind);
	free_and_null((char **) &qmatval);

	return(status);
}

static int set_quad_problem(double *H, double* P, double *C, double *B, char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p,char **sense_p, int **matbeg_p, int **matcnt_p,int **matind_p, double **matval_p, double **lb_p,double **ub_p, int **qmatbeg_p, int **qmatcnt_p, int **qmatind_p, double **qmatval_p){

	char     *zprobname = NULL;
	double   *zobj = NULL;
	double   *zrhs = NULL;
	char     *zsense = NULL;
	int      *zmatbeg = NULL;
	int      *zmatcnt = NULL;
	int      *zmatind = NULL;
	double   *zmatval = NULL;
	double   *zlb = NULL;
	double   *zub = NULL;
	int      *zqmatbeg = NULL;
	int      *zqmatcnt = NULL;
	int      *zqmatind = NULL;
	double   *zqmatval = NULL;
	int      status = 0;

	int i, j, k, l, m, ind, no_zeros = 0, no_qzeros = 0;
	double eps_val = 1e-9;

	//Computing the amount of non zeros elements in C matrix.
	for(i = 0;i < nc;i++){
		for(j = 0;j < 6*N;j++){
			if (fabs(*(C + 6*N*i + j)) > eps_val)
				no_zeros += 1;
		}
	}

	//Computing the amount of non zeros elements in H matrix.
	for(i = 0; i < N; i++){
		no_qzeros += 36*solid[i].cantBlock;
	}

	zprobname = (char *) malloc (16*sizeof(char));
	zobj      = (double *) malloc (6*N*sizeof(double));
	zrhs      = (double *) malloc (nc*sizeof(double));
	zsense    = (char *) malloc (nc*sizeof(char));
	zmatbeg   = (int *) malloc (6*N*sizeof(int));
	zmatcnt   = (int *) malloc (6*N*sizeof(int));
	zmatind   = (int *) malloc (no_zeros*sizeof(int));
	zmatval   = (double *) malloc (no_zeros*sizeof(double));
	zlb       = (double *) malloc (6*N*sizeof(double));
	zub       = (double *) malloc (6*N*sizeof(double));
	zqmatbeg  = (int *) malloc (6*N*sizeof(int));
	zqmatcnt  = (int *) malloc (6*N*sizeof(int));
	zqmatind  = (int *) malloc (no_qzeros * sizeof(int));
	zqmatval  = (double *) malloc (no_qzeros * sizeof(double));

	if ( zprobname == NULL || zobj    == NULL || zrhs      == NULL || zsense  == NULL || zmatbeg   == NULL || zmatcnt == NULL || zmatind   == NULL || zmatval == NULL || zlb == NULL || zub == NULL || zqmatbeg  == NULL ||zqmatcnt == NULL || zqmatind  == NULL || zqmatval == NULL)  {
		status = 1;
		goto TERMINATE;
	}

	strcpy(zprobname, "minimizing");

	//Building quadratic part
	m = 0;	
	for(i = 0; i < N; i++){
		for(j = 0; j < 6; j++){
			zqmatbeg[6*i + j] = m;
			zqmatcnt[6*i + j] = 6*solid[i].cantBlock;
			for(k = 0; k < solid[i].cantBlock; k++){
				ind = 36*N*solid[i].idBlock[k] + 6*i;
				for(l = 0; l < 6; l++){
					zqmatval[m] = H[ind + 6*N*l + j];
					zqmatind[m] = 6*solid[i].idBlock[k] + l;
					m++;
				}
			}			
		}
	}

	//Building linear part and constraint vectors
	k = 0;
	for(i = 0; i < 6*N; i++){
		zobj[i] = -P[i];
		zlb[i] = -CPX_INFBOUND;
		zub[i] = CPX_INFBOUND;
		l = 0;
		for(j = 0; j < nc; j++){
			if (fabs(C[6*N*j + i]) > eps_val){
				zmatval[k] = -C[6*N*j + i];
				zmatind[k] = j;
				k++;
				l++;
			}
		}
		if (i == 0)
			zmatbeg[0] = 0;
		else
			zmatbeg[i] = zmatbeg[i - 1] + zmatcnt[i - 1];
		zmatcnt[i] = l;
	}

	//Constraints right side --> -C.x <= -B
	for(i = 0; i < nc; i++){
		zsense[i] = 'L';
		zrhs[i] = -B[i];
	}

TERMINATE:
	if ( status ) {
		free_and_null ((char **) &zprobname);
		free_and_null ((char **) &zobj);
		free_and_null ((char **) &zrhs);
		free_and_null ((char **) &zsense);
		free_and_null ((char **) &zmatbeg);
		free_and_null ((char **) &zmatcnt);
		free_and_null ((char **) &zmatind);
		free_and_null ((char **) &zmatval);
		free_and_null ((char **) &zlb);
		free_and_null ((char **) &zub);
		free_and_null ((char **) &zqmatbeg);
		free_and_null ((char **) &zqmatcnt);
		free_and_null ((char **) &zqmatind);
		free_and_null ((char **) &zqmatval);
	}
	else {
		*numcols_p   = 6*N;
		*numrows_p   = nc;
		*objsen_p    = CPX_MIN;   //The problem is minimization

		*probname_p  = zprobname;
		*obj_p       = zobj;
		*rhs_p       = zrhs;
		*sense_p     = zsense;
		*matbeg_p    = zmatbeg;
		*matcnt_p    = zmatcnt;
		*matind_p    = zmatind;
		*matval_p    = zmatval;
		*lb_p        = zlb;
		*ub_p        = zub;
		*qmatbeg_p   = zqmatbeg;
		*qmatcnt_p   = zqmatcnt;
		*qmatind_p   = zqmatind;
		*qmatval_p   = zqmatval;
   }
   return (status);
}
