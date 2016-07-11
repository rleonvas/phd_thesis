#include "functions.h"
#include "vector.h"
#include "files.h"
#include "set_geometry.h"
#include "computing.h"

int main(int argc, char *argv[]){

	int i, it, k;
	char name_evol[20] = "evolution_", name_time[20] = "time_", ext[6] = ".data";
	contact *current;
  clock_t start,end;
  double cpu_time_used;
  FILE *file_t;
	
	strcat(strcat(name_evol,argv[1]),ext);
	file_evol = fopen(name_evol,"w");

	strcat(strcat(name_time,argv[1]),ext);
	file_t = fopen(name_time,"w");
	
	for (k = 0; k < 1; k++){
		start = clock();
		get_size();
		get_position();
		set_edges();
		set_faces();
		set_edge_faces();
		computing_vertices();
		computing_faces();
		computing_mass();	
		get_parameters();
		computing_inertia();
		computing_radius();

		computing_inertia_matrix();
		initialValues();
		it = 0;
		while (it < N_ITER){        
		  store_data();
			printf("TIMESTEP %d\n",it);
			computing_points();		
			if (it%ITPROX == 0){
		    	computing_proximity();					
			}
			/*Collision detection*/	
			current = contact_detection();

			if (current != NULL){				
        minimization(current);
			}		
			else{
				updating_velocity();
		  }
			updating_position();
			updating_axis();		
			computing_inertia_matrix();
			it++;		
		}       
		end = clock();
		cpu_time_used = ((double) (end - start))/CLOCKS_PER_SEC;        
		fprintf(file_t,"%.5f\n",cpu_time_used);        
	}
    fclose(file_t);
    fclose(file_evol);		
	return 0;
}
