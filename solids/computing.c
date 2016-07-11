#include "functions.h"
#include "set_geometry.h"
#include "computing.h"
#include "vector.h"

void computing_vertices(){

	int i = 0,j,k;
	double a = 1,b = -1,c = -1;

	while(i < 8){
		coef_ver[i][0] = a;
		coef_ver[i][1] = b;
		coef_ver[i][2] = c;
		if(i%4 == 3) a = (-1)*a;
		if(i%2 == 0) b = (-1)*b;
		if(i%2 == 1) c = (-1)*c;
		i++;
	}
	if(NVERT == 10){
		coef_ver[8][0] = -2.;coef_ver[8][1] = 0.;coef_ver[8][2] = 0.;
		coef_ver[9][0] = 2.;coef_ver[9][1] = 0.;coef_ver[9][2] = 0.;
	}
	for(i = 0; i < N; i++){
		for(j = 0; j < NVERT; j++){
			for(k = 0; k < 3; k++)
				solid[i].vertices[j][k] = 0.5*coef_ver[j][k]*solid[i].size[k];
		}
	}
}

void computing_points(){

	int i,j,k,l;
	double factor;

	for(i = 0; i < N; i++){
		for(j = 0; j < NVERT; j++){
			for(k = 0; k < 3; k++)
				solid[i].points[j][k] = solid[i].pos[k] + solid[i].vertices[j][0]*solid[i].axis[0][k] + solid[i].vertices[j][1]*solid[i].axis[1][k] + solid[i].vertices[j][2]*solid[i].axis[2][k];
		}
  	/*if (NDELTA > 0){
		  for(j = 0; j < NEDGE; j++){
		    for(k = 0; k < NDELTA; k++){
		      factor = (double)(k + 1)/(NDELTA + 1);
		      for(l = 0; l < 3; l++)
		        solid[i].points[NVERT + NDELTA*j + k][l] = solid[i].points[edge[j][0]][l] + factor*(solid[i].points[edge[j][1]][l] - solid[i].points[edge[j][0]][l]);
		    }
		  }
		}*/
	}
}

void computing_faces(){

	int i,j,k,l;
	double vec_one[3],vec_two[3],vec_three[3],vec_cross[3];
	double value;

	for(i = 0; i < N; i++){
		for(j = 0; j < NFACE; j++){
			for(k = 0; k < 3; k++){
				vec_one[k] = solid[i].vertices[face[j][1]][k] - solid[i].vertices[face[j][0]][k];
				vec_two[k] = solid[i].vertices[face[j][2]][k] - solid[i].vertices[face[j][0]][k];
				vec_three[k] = solid[i].vertices[face[j][3]][k];
			}
			cross_product(vec_one,vec_two,vec_cross);
			value = dot_product(vec_cross,vec_cross,3);
			for(k = 0; k < 3; k++){
				vec_cross[k] = vec_cross[k]/sqrt(value);
				solid[i].faces[j][k] = vec_cross[k];
			}
			value = dot_product(vec_cross,vec_three,3);
			solid[i].faces[j][3] = -value;
			for(l = 0; l < 4; l++){
			//printf("\nFace %d - Value %d --> %f",j,l,solid[i].faces[j][l]);
			}
		}		
	}
}

void computing_mass(){

	int i;

	for(i = 0; i < N; i++)
		solid[i].mass = DENSITY*solid[i].size[0]*solid[i].size[1]*solid[i].size[2];
}

void computing_inertia(){

	int i;
	double den = 12.;

	for(i = 0; i < N; i++){
		solid[i].diag_inertia[0] = solid[i].mass*(solid[i].size[1]*solid[i].size[1] + solid[i].size[2]*solid[i].size[2])/den;
		solid[i].diag_inertia[1] = solid[i].mass*(solid[i].size[0]*solid[i].size[0] + solid[i].size[2]*solid[i].size[2])/den;
		solid[i].diag_inertia[2] = solid[i].mass*(solid[i].size[0]*solid[i].size[0] + solid[i].size[1]*solid[i].size[1])/den;
	}
}

void computing_inertia_matrix(){

	int i,j,k,l;
	double aux[3][3];

	for(i = 0; i < N ; i++){
		for(j = 0; j < 3 ; j++){
			aux[0][j] = solid[i].axis[j][0]*solid[i].diag_inertia[0];
			aux[1][j] = solid[i].axis[j][1]*solid[i].diag_inertia[1];
			aux[2][j] = solid[i].axis[j][2]*solid[i].diag_inertia[2];
		}

    	for(j = 0; j < 3 ; j++){
	    	for(k = 0; k < 3 ; k++){
				solid[i].inertia_matrix[k][j] = 0.0;
				for(l = 0;l < 3; l++)
solid[i].inertia_matrix[k][j] = solid[i].inertia_matrix[k][j] + aux[l][j]*solid[i].axis[k][l];
			}
		}
	}
}

void computing_radius(){

	int i,j;
	double r;

	for(i = 0; i < N; i++){
		r = 0.;
		for(j = 0;j < 3; j++)
			r += solid[i].size[j]*solid[i].size[j];
		solid[i].rad = sqrt(r);
	}
}

