#include "functions.h"
#include "vector.h"

void initialValues(){
	int i;

	for (i = 0; i < N; i++){
		solid[i].C1 = (solid[i].diag_inertia[1] - solid[i].diag_inertia[2])/solid[i].diag_inertia[0];
		solid[i].C2 = (solid[i].diag_inertia[2] - solid[i].diag_inertia[0])/solid[i].diag_inertia[1];
		solid[i].C3 = (solid[i].diag_inertia[0] - solid[i].diag_inertia[1])/solid[i].diag_inertia[2];
	}
}

void free_and_null(char **ptr){
	if (*ptr != NULL){
		free(*ptr);
		*ptr = NULL;
	}
}

void displayMatrix(double *mat, int dim_row, int dim_col){
	int i,j;

	for(i = 0; i < dim_row; i++){
		for(j = 0; j < dim_col; j++)
			printf("%.9f, ",mat[dim_col*i + j]);
		printf("\n");
	}
	printf("\n");
}

contact* contact_detection(){
	int lf, lc;
  contact *contact_pt = NULL;

  nc = 0;

	for(lf = 0; lf < N; lf++){
		solid[lf].cantBlock = 0;
		for(lc = 0; lc < 30; lc++)
			solid[lf].idBlock[lc] = -1;
	}

  contact_pt = surface_contact_detection(contact_pt);
  contact_pt = solid_contact_detection(contact_pt);

  return contact_pt;
}

contact* new_contact(contact* list, int sol_one, int sol_two, double gp_one[3], double gp_two[3], double normal_vec[3]){
	int i;
	contact* current;

	current = malloc(sizeof(contact));
	current->sol_one = sol_one;
	current->sol_two = sol_two;
	for(i = 0; i < 3; i++){
		current->gp_one[i] = gp_one[i];
		current->gp_two[i] = gp_two[i];
		current->normal_vec[i] = normal_vec[i];
	}
	current->next = list;
	list = current;
	nc += 1;
	
  return list;
}

int sign(double x){    
    if (x >= 0){
        return 1;
    }
    else{
        return -1;
    }
}

int grid_contact(double point[3], double normal_vec[3]){
	int ind = 0, i ;
	double norm;
	
	if (point[2] <= 0){
		ind = 1;
		normal_vec[0] = 0.;
		normal_vec[1] = 0.;
		normal_vec[2] = 1.;
	}
	return ind;
}


contact* surface_contact_detection(contact *list){
	int i,j,k,ind;
	double test_point[3],normal_vec[3],vec_cross[3];
	double ga_point[3],rel_vel[3], gp_zeros[] = {0.,0.,0.};
	double test_vel;

	for(i = 0; i < N; i++){		
		for(j = 0; j < NVERT; j++){
			for(k = 0;k < 3;k++)
				test_point[k] = solid[i].points[j][k];
			//ind = 1 --> Vertex is below the surface.
			ind = grid_contact(test_point,normal_vec);
			if (ind == 1){
				//ga_point --> Contact point to solid position.
				for(k = 0; k < 3; k++)
					ga_point[k] = test_point[k] - solid[i].pos[k];
				//Cross product --> solid.omg x ga_point
				cross_product(solid[i].omg,ga_point,vec_cross);
				//Relative velocity.
				for(k = 0; k < 3; k++)
					rel_vel[k] = solid[i].vel[k] + vec_cross[k];
				test_vel = dot_product(rel_vel,normal_vec,3);
				//Testing if the solid moves inwards or outwards from the surface.
				if (test_vel <= 0)
					list = new_contact(list,i,-1,ga_point,gp_zeros,normal_vec);
			}				
		}		
	}
	return list;
}

int overlapping(int sol_i, int sol_j){
	int lf;
	double rad_i = 0., rad_j = 0., d_pos = 0.;

	for(lf = 0;lf < 3; lf++)
		d_pos +=(solid[sol_i].pos[lf]-solid[sol_j].pos[lf])*(solid[sol_i].pos[lf]-solid[sol_j].pos[lf]);
	d_pos = sqrt(d_pos);

	if (d_pos > (solid[sol_i].rad + solid[sol_j].rad))
		return 0;
	else
		return 1;
}

int point_inside_polygon(double point[3], int sol_ind){
	/*Check if point point_p is inside of solid sol_ind.*/
	int i, j, status = 1;
	double point_pos[3],point_ref[3], dist, eps_value = 1e-12;

	/*Point must be in the sol_ind solid coordinates*/
	for(i = 0;i < 3; i++)
		point_pos[i] = point[i] - solid[sol_ind].pos[i];
	for(i = 0;i < 3; i++){
		point_ref[i] = 0.;
		for(j = 0; j < 3; j++)
			point_ref[i] += point_pos[j]*solid[sol_ind].axis[i][j];
  }
	/*Computing distance to each polygon face, if the point is inside, then
	all distances must be negative. There is a break if any distance is positive.*/	
	for(i = 0; i < NFACE; i++){				
		dist = solid[sol_ind].faces[i][3];
		for(j = 0; j < 3; j++)
			dist += point_ref[j]*solid[sol_ind].faces[i][j];
		if (dist > eps_value){
			status = 0;
			break;
		}
	}	
	return status;
}

int face_min_distance(double point[3], int sol_point, int sol_face){
	/*Return the face which is closer to the point point_p.*/
	int i, j, ind;
	double point_con[3], point_pos[3], point_ref_con[3], point_ref_sol[3];
  double dist, test_value, min_dist = 1e+12, eps_value = 1e-12;

	/*Point must be in solid sol_ind coordinates*/
	for(i = 0;i < 3; i++){
		point_con[i] = point[i] - solid[sol_face].pos[i];
		point_pos[i] = solid[sol_point].pos[i] - solid[sol_face].pos[i];
	}
	for(i = 0;i < 3; i++){
		point_ref_con[i] = 0.;
    point_ref_sol[i] = 0.;
		for(j = 0; j < 3; j++){
			point_ref_con[i] += point_con[j]*solid[sol_face].axis[i][j];
      point_ref_sol[i] += point_pos[j]*solid[sol_face].axis[i][j];
		}
	}

	/*Computing each distance face and obtaining the minimal one.*/
	for(i = 0; i < NFACE; i++){
		dist = solid[sol_face].faces[i][3];
    test_value = 0.;
		for(j = 0; j < 3; j++){
		  dist += point_ref_con[j]*solid[sol_face].faces[i][j];
		  test_value += (point_ref_con[j] - point_ref_sol[j])*solid[sol_face].faces[i][j];
		}
    if((fabs(dist) < fabs(min_dist)) && (test_value < 0.)){
      ind = i;
      min_dist = dist;
		}
	}
	return ind;
}

contact* solid_contact_detection(contact *list){

	int i, j, k, l, m, ind, status, lf, lc, t, ver_ini, ver_end, coEdge, face_one, face_two;
	double test_vel, norm, eps_value = 1e-12, rel_vel[3], normal_vec[3];
	/*----------------- face contact detection variables -----------------*/
	double test_point[3], ref_point[3], ref_face[3];
	double gp_point[3], gp_face[3];	
	/*----------------- edge contact detection variables -----------------*/
	double vec_edge_ini[3], vec_edge_end[3], vec_ini[3], vec_end[3], ref_ini[3], ref_end[3], c_point[3];
	double s_ini[3], s_end[3], s_min, s_max, js[3][2], s[2], s_ave, diff_ref;
	double gp_edge_i[3], gp_edge_j[3], cp_one[3],cp_two[3], cEdge_one[3],cEdge_two[3];
	double rel_omg_i[3], rel_omg_j[3], vec_faces[3], sum_faces[3], edgeNormal;

	for(i = 0; i < N; i++){
		if (solid[i].idBlock[solid[i].cantBlock] != -1) solid[i].cantBlock++;
		solid[i].idBlock[solid[i].cantBlock] = i;
		for(t = 0; t < solid[i].n_test; t++){
      j = solid[i].tabprox[t];
			if ((i < j) && (overlapping(i,j))){
				/*--------------FACE CONTACT DETECTION-------------------*/
        for(k = 0; k < NPOINTS; k++){
          for(l = 0;l < 3; l++)
          	test_point[l] = solid[i].points[k][l];
					status = point_inside_polygon(test_point,j);
					if(status){
			      ind = face_min_distance(test_point,i,j);
						for(l = 0; l < 3; l++){
				      gp_point[l] = test_point[l] - solid[i].pos[l];
							gp_face[l] = test_point[l] - solid[j].pos[l];
							normal_vec[l] = 0.;
							for(m = 0; m < 3; m++)
			          normal_vec[l] += solid[j].faces[ind][m]*solid[j].axis[m][l];
			    	}
				    norm = sqrt(dot_product(normal_vec,normal_vec,3));
						cross_product(solid[i].omg,gp_point,ref_point);
						cross_product(solid[j].omg,gp_face,ref_face);
						for(l = 0; l < 3; l++){
				      normal_vec[l] /= norm;
				      rel_vel[l] = solid[i].vel[l] + ref_point[l] - solid[j].vel[l] - ref_face[l];
				    }
						test_vel = dot_product(rel_vel,normal_vec,3);
						if (test_vel <= eps_value){          
							list = new_contact(list,i,j,gp_point,gp_face,normal_vec);
							if (solid[i].idBlock[solid[i].cantBlock] != j){
								solid[i].cantBlock++;
								solid[i].idBlock[solid[i].cantBlock] = j;
							}
							for(m = 0; m < 3; m++)
							  normal_vec[m] = -1.*normal_vec[m];
							list = new_contact(list,j,i,gp_face,gp_point,normal_vec);
							if (solid[j].idBlock[solid[j].cantBlock] != i){								
								if(solid[j].idBlock[solid[j].cantBlock] != -1) solid[j].cantBlock++;
								solid[j].idBlock[solid[j].cantBlock] = i;
							}
							//printf("CONTACT %d - SOLID %d,SOLID %d\n",nc,i,j);
						}
					}
				}
				/*--------------FACE CONTACT DETECTION-------------------*/

				/*--------------EDGE CONTACT DETECTION-------------------*/
				for(k = 0; k < NEDGE; k++){
					status = 1;
					ver_ini = edge[k][0];
					ver_end = edge[k][1];
					for(lf = 0; lf < 3; lf++){
						vec_edge_ini[lf] = solid[j].points[ver_ini][lf];
						vec_edge_end[lf] = solid[j].points[ver_end][lf];
					}
					//Edge k - Solid j in reference system of Solid i
					for(lf = 0; lf < 3; lf++){
						vec_ini[lf] = vec_edge_ini[lf] - solid[i].points[4][lf];
						vec_end[lf] = vec_edge_end[lf] - solid[i].points[4][lf];
					}
					for(lf = 0; lf < 3; lf++){
						ref_ini[lf] = 0;
						ref_end[lf] = 0;
						for(lc = 0; lc < 3; lc++){
							ref_ini[lf] += solid[i].axis[lf][lc]*vec_ini[lc];
							ref_end[lf] += solid[i].axis[lf][lc]*vec_end[lc];
						}
					}
					//Computing parameter s_min and s_max for each position
					for(lf = 0; lf < 3; lf++){
						diff_ref = (ref_end[lf] - ref_ini[lf]);
						if (fabs(diff_ref) < eps_value){
							//Checking if it is in the interval [0,size_lf]
							if (ref_ini[lf] >= 0 && ref_ini[lf] <= solid[i].size[lf]){
								s_min = 0.;
								s_max = 1.;
							}
							else{
								status = 0;
								break;
							}
						}
						else{
							s_ini[lf] = -ref_ini[lf]/diff_ref;
							s_end[lf] = (solid[i].size[lf] - ref_ini[lf])/diff_ref;
							s_min = (s_ini[lf] < s_end[lf]) ? s_ini[lf] : s_end[lf];
							s_max = s_ini[lf] + s_end[lf] - s_min;
						}
						js[lf][0] = (s_min > 0.) ? s_min : 0.;
						js[lf][1] = (s_max < 1.) ? s_max : 1.;
						if (js[lf][0] > js[lf][1]){ 
							status = 0;
							break;
						}
					}
					if (status){					
						s[0] = max_value(js[0][0],js[1][0],js[2][0]);
						s[1] = min_value(js[0][1],js[1][1],js[2][1]);
						if (s[0] < s[1]){
							s_ave = .5*(s[0] + s[1]);
							for (lf = 0; lf < 3; lf++){
								cp_one[lf] = (1 - s[0])*vec_edge_ini[lf] + s[0]*vec_edge_end[lf];
								cp_two[lf] = (1 - s[1])*vec_edge_ini[lf] + s[1]*vec_edge_end[lf];
								c_point[lf] = (1 - s_ave)*vec_edge_ini[lf] + s_ave*vec_edge_end[lf];
							}
							//Obtain the contact edge of solid i.
							coEdge = getContactEdge(cp_one,cp_two,i);						
							face_one = edge_faces[coEdge][0];
							face_two = edge_faces[coEdge][1];
							//Computing face_one + face_two
							for(lf = 0; lf < 3; lf++){
								sum_faces[lf] = solid[i].faces[face_one][lf] + solid[i].faces[face_two][lf];
							}
							for(lf = 0; lf < 3; lf++){
								vec_faces[lf] = 0;
								for(lc = 0; lc < 3; lc++){
									vec_faces[lf] += solid[i].axis[lc][lf]*sum_faces[lc];
								}
							}

							for (lf = 0; lf < 3; lf++){
								gp_edge_i[lf] = c_point[lf] - solid[i].pos[lf];
								gp_edge_j[lf] = c_point[lf] - solid[j].pos[lf];
								cEdge_one[lf] = solid[i].points[edge[coEdge][1]][lf] - solid[i].points[edge[coEdge][0]][lf];
								cEdge_two[lf] = vec_edge_end[lf] - vec_edge_ini[lf];
							}
						
							//Computing normal vector edge_i x edge_j
							cross_product(cEdge_one, cEdge_two, normal_vec);
							norm = sqrt(dot_product(normal_vec,normal_vec,3));
							for(lf = 0; lf < 3; lf++)
								normal_vec[lf] /= norm;
							edgeNormal = dot_product(vec_faces,normal_vec,3);
							if (edgeNormal > 0){
								for(lf = 0; lf < 3; lf++)
									normal_vec[lf] = -normal_vec[lf];
							}
							//Computing relative rotational velocity
							cross_product(solid[i].omg,gp_edge_i,rel_omg_i);
							cross_product(solid[j].omg,gp_edge_j,rel_omg_j);
							//Computing total relative velocity
							for(lf = 0; lf < 3; lf++)
								rel_vel[lf] = solid[i].vel[lf] + rel_omg_i[lf] - solid[j].vel[lf] - rel_omg_j[lf];
							//Verifying inwards or outwards velocity
							test_vel = dot_product(rel_vel,normal_vec,3);
							if (test_vel <= eps_value){
								list = new_contact(list,i,j,gp_edge_i,gp_edge_j,normal_vec);
								if (solid[i].idBlock[solid[i].cantBlock] != j){								
									solid[i].cantBlock++;
									solid[i].idBlock[solid[i].cantBlock] = j;
								}
								for (lf = 0; lf < 3; lf++)
										normal_vec[lf] = -normal_vec[lf];
								list = new_contact(list,j,i,gp_edge_j,gp_edge_i,normal_vec);
								if (solid[j].idBlock[solid[j].cantBlock] != i){								
									if(solid[j].idBlock[solid[j].cantBlock] != -1) solid[j].cantBlock++;
									solid[j].idBlock[solid[j].cantBlock] = i;
								}
								//printf("CONTACT %d - SOLID %d,SOLID %d\n",nc,i,j);
							}
						}
					}
				}
				/*--------------EDGE CONTACT DETECTION-------------------*/
			}
		}
		solid[i].cantBlock++;
	}
	return list;
}

int getContactEdge(double cp_one[3], double cp_two[3], int idSolid){
	int i,j,face_one,face_two;
	double evalFace, ref_one[3], ref_two[3];
	//Point one and two in local reference system of solid i.
	for(i = 0; i < 3; i++){
		cp_one[i] -= solid[idSolid].pos[i];
		cp_two[i] -= solid[idSolid].pos[i];
	}
	for(i = 0; i < 3; i++){
		ref_one[i] = 0.;
		ref_two[i] = 0.;
		for(j = 0; j < 3; j++){
			ref_one[i] += solid[idSolid].axis[i][j]*cp_one[j];
			ref_two[i] += solid[idSolid].axis[i][j]*cp_two[j];
		}
	}
	//Computing the contact faces.
	for(i = 0; i < NFACE; i++){
		evalFace = 0.;
		for(j = 0; j < 3; j++)
			evalFace += solid[idSolid].faces[i][j]*ref_one[j];		
		evalFace += solid[idSolid].faces[i][3];
		if (fabs(evalFace) <= 1e-12){
			face_one = i;
			break;
		}
	}
	for(i = 0; i < NFACE; i++){
		evalFace = 0.;
		for(j = 0; j < 3; j++)
			evalFace += solid[idSolid].faces[i][j]*ref_two[j];		
		evalFace += solid[idSolid].faces[i][3];
		if (fabs(evalFace) <= 1e-12){
			face_two = i;
			break;
		}
	}	
	//Searching edge associated to face_one/face_two
	for(i = 0; i < NEDGE; i++){
		if ((edge_faces[i][0] == face_one && edge_faces[i][1] == face_two) || (edge_faces[i][0] == face_two && edge_faces[i][1] == face_one))
			return i;
	}	
}

void updating_velocity(){
	int i,j;	
	double omg1,omg2,omg3;
	
	for (i = 0;i < N;i++){
		//Updating velocities without contacts.
		solid[i].vel[2] = solid[i].vel[2] - G*timestep;
		//Updating rotational velocity with Euler.
		omg1 = solid[i].omg[0];
		omg2 = solid[i].omg[1];
		omg3 = solid[i].omg[2];
		solid[i].omg[0] = omg1 + timestep*solid[i].C1*omg2*omg3;
		solid[i].omg[1] = omg2 + timestep*solid[i].C2*omg1*omg3;
		solid[i].omg[2] = omg3 + timestep*solid[i].C3*omg1*omg2;
	}
}

void updating_position(){
	int i,j;

	for(i = 0;i < N;i++){
		solid[i].pos[0] += solid[i].vel[0]*timestep;
		solid[i].pos[1] += solid[i].vel[1]*timestep;
		solid[i].pos[2] += solid[i].vel[2]*timestep - .5*G*timestep*timestep;
	}
}

void updating_axis(){

	int i,j,k,l;
	double *ROT,ia[3],ja[3],ka[3],angle_rot[3];
	double norm_ia,norm_ja,norm_ka;

	ROT = (double*)calloc(9,sizeof(double));
	for(i = 0; i < N; i++){
		//Rotation angle.
		angle_rot[0] = solid[i].omg[0]*timestep;
		angle_rot[1] = solid[i].omg[1]*timestep;
		angle_rot[2] = solid[i].omg[2]*timestep;

		//Computing rotation matrix.		
		*(ROT) = cos(angle_rot[1])*cos(angle_rot[2]);
		*(ROT + 1) = sin(angle_rot[0])*sin(angle_rot[1])*cos(angle_rot[2]) - cos(angle_rot[0])*sin(angle_rot[2]);
		*(ROT + 2) = cos(angle_rot[0])*sin(angle_rot[1])*cos(angle_rot[2]) + sin(angle_rot[0])*sin(angle_rot[2]);
		*(ROT + 3) = cos(angle_rot[1])*sin(angle_rot[2]);
		*(ROT + 4) = sin(angle_rot[0])*sin(angle_rot[1])*sin(angle_rot[2]) + cos(angle_rot[0])*cos(angle_rot[2]);
		*(ROT + 5) = cos(angle_rot[0])*sin(angle_rot[1])*sin(angle_rot[2]) - sin(angle_rot[0])*cos(angle_rot[2]);
		*(ROT + 6) = -sin(angle_rot[1]);
		*(ROT + 7) = sin(angle_rot[0])*cos(angle_rot[1]);
		*(ROT + 8) = cos(angle_rot[0])*cos(angle_rot[1]);

		//Updating axis of the solid.
		for(j = 0; j < 3; j++){
			ia[j] = 0.;
			ja[j] = 0.;
			ka[j] = 0.;
			for(k = 0; k < 3; k++){
				ia[j] += ROT[3*j + k]*solid[i].axis[0][k];
				ja[j] += ROT[3*j + k]*solid[i].axis[1][k];
				ka[j] += ROT[3*j + k]*solid[i].axis[2][k];
			}
		}
		norm_ia = dot_product(ia, ia, 3);
		norm_ja = dot_product(ja, ja, 3);
		norm_ka = dot_product(ka, ka, 3);
		for(j = 0; j < 3;j++){
			ia[j] = ia[j]/sqrt(norm_ia);
			ja[j] = ja[j]/sqrt(norm_ja);
			ka[j] = ka[j]/sqrt(norm_ka);
		}
		for(j = 0; j < 3;j++){
			solid[i].axis[0][j] = ia[j];
			solid[i].axis[1][j] = ja[j];
			solid[i].axis[2][j] = ka[j];
		}		
	}
	free_and_null((char **) &ROT);
}

void store_data(){
  int i,j,k;

	for(i = 0; i < N ; i++){
		for(j = 0; j < 3; j++)
			fprintf(file_evol,"%1.9f ",solid[i].pos[j]);
		fprintf(file_evol,"\n");
		for(j = 0; j < 3; j++)
			fprintf(file_evol,"%1.9f ",solid[i].vel[j]);
		fprintf(file_evol,"\n");
		for(j = 0; j < 3; j++)
			fprintf(file_evol,"%1.9f ",solid[i].omg[j]);
		fprintf(file_evol,"\n");
		for(j = 0; j < 3; j++){
			for(k = 0; k < 3; k++)
				fprintf(file_evol,"%1.9f ",solid[i].axis[j][k]);	
			fprintf(file_evol,"\n");
		}
	}
}

void computing_proximity(){
  int i, j, k, l, m, n, ii, jj, n_test, id_solid, px, py, pz, di, dj, dk, nb, *listNB = NULL, max_nb;
	int diff_pos[26][3] = {{-1,-1,-1,},{-1,-1,0},{-1,-1,1},{-1,0,-1},{-1,0,0},{-1,0,1},{-1,1,-1},{-1,1,0},{-1,1,1},{0,-1,-1},{0,-1,0},{0,-1,1},{0,0,-1},{0,0,1},{0,1,-1},{0,1,0},{0,1,1},{1,-1,-1},{1,-1,0},{1,-1,1},{1,0,-1},{1,0,0},{1,0,1},{1,1,-1},{1,1,0},{1,1,1}};
  double eps_pos = 1e-6, coef_d = 2.*MAX_S, pfx, pfy, pfz;
	nbhood* array;
	double* lim_pos;

	lim_pos = position_limit();

	ddx = (int)(ceil((eps_pos + lim_pos[1] - lim_pos[0] - eps_pos)/coef_d));
	ddy = (int)(ceil((eps_pos + lim_pos[3] - lim_pos[2] - eps_pos)/coef_d));
	ddz = (int)(ceil((eps_pos + lim_pos[5] - lim_pos[4] - eps_pos)/coef_d));

	/*printf("ddx = %d ",ddx);
	printf("ddy = %d ",ddy);
	printf("ddz = %d \n",ddz);*/

	array = (nbhood*)malloc(ddx*ddy*ddz*sizeof(nbhood));
	for(i = 0; i < ddx; i++){
		for(j = 0; j < ddy; j++){
			for(k = 0; k < ddz; k++)
				DD(i,j,k).count = 0;
		}
	}

	for(i = 0; i < N; i++){
		px = (int)(floor((solid[i].pos[0] - lim_pos[0])/coef_d));
		py = (int)(floor((solid[i].pos[1] - lim_pos[2])/coef_d));
		pz = (int)(floor((solid[i].pos[2] - lim_pos[4])/coef_d));
		DD(px,py,pz).l_of_sol[DD(px,py,pz).count] = i;
		DD(px,py,pz).count++;
	}

	//Computing maximum number of neighboors.
	max_nb = 0;
	for(i = 0; i < ddx; i++){
		for(j = 0; j < ddy; j++){
			for(k = 0; k < ddz; k++){
				/*if (DD(i,j,k).count > 0){
					printf("(%d,%d,%d)(%d): ",i,j,k,ddy*ddz*i + ddz*j + k);
					for(l = 0; l < DD(i,j,k).count; l++)
						printf("%d ",DD(i,j,k).l_of_sol[l]);
					printf("\n");
				}*/
				if(DD(i,j,k).count > max_nb) max_nb = DD(i,j,k).count;
			}				
		}
	}
	//printf("MAX NB: %d\n\n",max_nb);	
	//Computing neighboor list for each solid	
	for(i = 0; i < ddx; i++){
		for(j = 0; j < ddy; j++){
			for(k = 0; k < ddz; k++){
				if (DD(i,j,k).count > 0){
					//printf("-------------<%d,%d,%d>-------------\n",i,j,k);
					//Computing the closer solids to the neighborhood DD(i,j,k)
					nb = 0;
					listNB = (int*)malloc(26*max_nb*sizeof(int));
					for(l = 0; l < 26; l++){						
						di = i + diff_pos[l][0];
						dj = j + diff_pos[l][1];
						dk = k + diff_pos[l][2];						
						if (inRange(di,0,ddx) && inRange(dj,0,ddy) && inRange(dk,0,ddz) && (DD(di,dj,dk).count>0)){
							//printf("l = %d <di,dj,dk> = <%d,%d,%d>\n",l,di,dj,dk);
							for(m = 0; m < DD(di,dj,dk).count; m++){
								//printf("%d ",DD(di,dj,dk).l_of_sol[m]);
								*(listNB + nb) = DD(di,dj,dk).l_of_sol[m];
								nb++;
							}
							//printf("\n");
						}						
					}
					/*printf("\nNB DD(%d,%d,%d): ",i,j,k);
					for(l = 0; l < nb; l++)
						printf("%d ",listNB[l]);*/
					//Computing for each solid in DD(i,j,k), the proximity neighborhood.
					for(ii = 0; ii < DD(i,j,k).count; ii++){
						//Adding the rest of DD(i,j,k) to tabprox solid.
						id_solid = DD(i,j,k).l_of_sol[ii];
						//printf("\nSOLID: %d",id_solid);
						n_test = 0;
						for(l = 0; l < DD(i,j,k).count; l++){
							if (DD(i,j,k).l_of_sol[l] != id_solid){
								solid[id_solid].tabprox[n_test] = DD(i,j,k).l_of_sol[l];
								n_test++;
							}
						}	
						//printf("\n");
						//Adding the closer solids to tabprox solid.					
						for(l = 0; l < nb; l++){
							//printf("list_nb(%d) = %d\n",l,listNB[l]);
							solid[id_solid].tabprox[n_test] = listNB[l];
							n_test++;
						}
						//printf("\nTABPROX\n");
						solid[id_solid].n_test = n_test;						
						merge_sort(solid[id_solid].tabprox,n_test);
						/*for(l = 0; l < n_test; l++)
							printf("%d ",solid[id_solid].tabprox[l]);
						printf("\n");*/
					}
					free(listNB);
				}				
			}				
		}
	}
	free(array);
}

double* position_limit(){
	int i;
	double xmin = 1e10, ymin = 1e10, zmin = 1e10;
	double xmax = -1e10, ymax = -1e10, zmax = -1e10;
	static double lim[6];

	for(i = 0; i < N; i++){
		if (solid[i].pos[0] < xmin) xmin = solid[i].pos[0];
		if (solid[i].pos[1] < ymin) ymin = solid[i].pos[1];
		if (solid[i].pos[2] < zmin) zmin = solid[i].pos[2];
		if (solid[i].pos[0] > xmax) xmax = solid[i].pos[0];
		if (solid[i].pos[1] > ymax) ymax = solid[i].pos[1];
		if (solid[i].pos[2] > zmax) zmax = solid[i].pos[2];
	}
	lim[0] = xmin;
	lim[1] = xmax;
	lim[2] = ymin;
	lim[3] = ymax;
	lim[4] = zmin;
	lim[5] = zmax;
	
	return lim;
}

int inRange(int x, int lb, int ub){
	return (x >= lb) && (x < ub);
}
