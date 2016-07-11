#include "functions.h"
#include "files.h"

void get_size(){
	double value;
	char name[20] = "size_", ext[10] = ".data", line[128], num[10], *tokptr = NULL;
	int i,j;
	FILE *fp;

	sprintf(num,"%d",N);
    strcat(strcat(name,num),ext);
	fp = fopen(name,"r");
	if (fp != NULL){
		i = 0;
		while(fgets(line,sizeof(line),fp) != NULL){
			j = 0;
			tokptr = strtok(line," ");
			while(tokptr != NULL){
				value = atof(tokptr);
				solid[i].size[j] = value;
				tokptr = strtok(NULL," ");
				j++;
			}
			//Set the id of the solid.
			solid[i].id = i;
			i++;
		}
		fclose(fp);
	}
	else{
		printf("Cannot open file\n");
	}
}

void get_position(){
	double value;
	char name[20] = "position_", ext[10] = ".data", line[128], num[10], *tokptr = NULL;
	int i,j,k,ind;
	FILE *fp;

	sprintf(num,"%d",N);
    strcat(strcat(name,num),ext);
	fp = fopen(name,"r");
	if (fp != NULL){
		i = 0;
		k = 0;
		while(fgets(line,sizeof(line),fp) != NULL){
			j = 0;
			tokptr = strtok(line," ");
			ind = i%6;
			while(tokptr != NULL){
				value = atof(tokptr);
				switch(ind){
					case 0:
						solid[k].pos[j] = value;
						break;
					case 1:
						solid[k].vel[j] = value;
						break;
					case 2:
						solid[k].omg[j] = value;
						break;
					case 3:
						solid[k].axis[0][j] = value;
						break;
					case 4:
						solid[k].axis[1][j] = value;
						break;
					case 5:
						solid[k].axis[2][j] = value;
						if(j == 2) k++;
						break;
				}
				tokptr = strtok(NULL," ");
				j++;
			}
			i++;
		}
		fclose(fp);
	}
	else{
		printf("Cannot open file\n");
	}
}

void get_parameters(){
    /*parameter.data timestep n_iter it_prox k_n k_t*/

	double total_mass = 0.0, value, coef_norm, coef_tang;
	int i;
	char line[128],*tokptr = NULL;
	FILE *fp;

	fp = fopen("parameter.data","r");
	if (fp != NULL){
		//Timestep.
		fgets(line,sizeof(line),fp);
		tokptr = strtok(line," ");
		value = atof(tokptr);
		timestep = value;
		tokptr = strtok(NULL," ");
		//Number of iterations.
		fgets(line,sizeof(line),fp);
		tokptr = strtok(line," ");
		value = atof(tokptr);
		N_ITER = value;
		tokptr = strtok(NULL," ");
		//Delta iteration for computing neighbors.
		fgets(line,sizeof(line),fp);
		tokptr = strtok(line," ");
		value = atof(tokptr);
		ITPROX = value;
		tokptr = strtok(NULL," ");
		//Normal parameter.
		fgets(line,sizeof(line),fp);
		tokptr = strtok(line," ");
		value = atof(tokptr);
		coef_norm = value;
		tokptr = strtok(NULL," ");
		//Tagential parameter.
		fgets(line,sizeof(line),fp);
		tokptr = strtok(line," ");
		value = atof(tokptr);
		coef_tang = value;
		tokptr = strtok(NULL," ");
		fclose(fp);
	}
	else{
		printf("Cannot open file\n");
	}	
	
	for(i = 0;i < N; i++){
		total_mass += solid[i].mass;
	}
	//printf("COEF NORM: %f\n",coef_norm);
	//printf("COEF TANG: %f\n",coef_tang);
	//Normal parameter
	kn = coef_norm*total_mass/(double)N;
	//Tangential parameter
	kt = coef_tang*total_mass/(double)N;	
}
