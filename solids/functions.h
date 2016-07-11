#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#define N 3000
#define G 9.80665
#define NDELTA 0
#define NVERT 8
#define NEDGE 12
#define NFACE 6
#define NPROX 2000
#define DENSITY 5000
#define MAX_S 0.5

#define NPOINTS NVERT + NDELTA*NEDGE

int ddx, ddy, ddz;
#define DD(i,j,k) (array[ddy*ddz*i + ddz*j + k])

typedef struct nbhodd{
	int l_of_sol[100];
	int count;
}nbhood;

typedef struct contact{
  int sol_one;		//ID solid 1
	int sol_two;		//ID solid 2
	double gp_one[3];	//Vector from center of solid 1 to contact point.
	double gp_two[3];	//Vector from center of solid 2 to contact point.
	double normal_vec[3];	//Normal vector in the contact point.
	struct contact *next;
}contact;

struct solid_structure{
	int id;				//identification number.
	double mass;			//mass of the solid.
	double size[3];			//size of the solid.
	double pos[3];			//current position.
	double vel[3];			//transation velocity.
	double omg[3];			//rotation velocity.
	double axis[3][3];		//axis i,j,k.
	double diag_inertia[3];		//diagonal inertia matrix.
	double inertia_matrix[3][3];	//inertia matrix.
	double vertices[NVERT][3];	//vertices.
	double points[NPOINTS][3];	//points in the general reference.
  double faces[NFACE][4];	        //planes equation for each face.
  int tabprox[NPROX];             //NPROX closest solids.
	int n_test;
	double C1,C2,C3;
	double rad;
	int idBlock[100];
	int cantBlock;
} solid[N];

double coef_ver[NVERT][3];
double kn,kt,timestep;
int edge[NEDGE][2],face[NFACE][4],edge_faces[NEDGE][2];
int nc,N_ITER,ITPROX;
FILE *file_evol;

void initialValues();

void free_and_null(char **ptr);

void displayMatrix(double *mat, int dim_row, int dim_col);

contact* contact_detection();

contact* new_contact(contact* list, int sol_one, int sol_two, double gp_one[3], double gp_two[3], double normal_vec[3]);

int sign(double x);

int grid_contact(double point[3], double normal_vec[3]);

contact* surface_contact_detection(contact *list);

int overlapping(int sol_i, int sol_j);

int point_inside_polygon(double point[3], int sol_ind);

int face_min_distance(double point[3], int sol_point, int sol_face);

contact* solid_contact_detection(contact *list);

int getContactEdge(double cp_one[3], double cp_two[3], int idSolid);

void updating_velocity();

void updating_position();

void updating_axis();

void store_data();

void computing_proximity();

double* position_limit();

int inRange(int x, int lb, int ub);

#endif
