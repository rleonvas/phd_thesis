#include "functions.h"
#include "set_geometry.h"

void set_edges(){

	edge[0][0] = 0;edge[0][1] = 1;edge[1][0] = 1;edge[1][1] = 2;
	edge[2][0] = 2;edge[2][1] = 3;edge[3][0] = 3;edge[3][1] = 0;
	edge[4][0] = 2;edge[4][1] = 6;edge[5][0] = 7;edge[5][1] = 3;
	edge[6][0] = 0;edge[6][1] = 4;edge[7][0] = 1;edge[7][1] = 5;
	edge[8][0] = 5;edge[8][1] = 4;edge[9][0] = 4;edge[9][1] = 7;
	edge[10][0] = 6;edge[10][1] = 7;edge[11][0] = 6;edge[11][1] = 5;
	if(NEDGE == 20){
		edge[12][0] = 0;edge[12][1] = 8;edge[13][0] = 3;edge[13][1] = 8;
		edge[14][0] = 4;edge[14][1] = 8;edge[15][0] = 7;edge[15][1] = 8;
		edge[16][0] = 1;edge[16][1] = 9;edge[17][0] = 2;edge[17][1] = 9;
		edge[18][0] = 5;edge[18][1] = 9;edge[19][0] = 6;edge[19][1] = 9;
	}
}

void set_faces(){
	face[0][0] = 0;face[0][1] = 1;face[0][2] = 2;face[0][3] = 3;
	face[1][0] = 3;face[1][1] = 2;face[1][2] = 6;face[1][3] = 7;
	face[2][0] = 7;face[2][1] = 6;face[2][2] = 5;face[2][3] = 4;
	face[3][0] = 4;face[3][1] = 5;face[3][2] = 1;face[3][3] = 0;
	if(NFACE == 12){
		face[4][0] = 0;face[4][1] = 8;face[4][2] = 4;face[4][3] = 0;
		face[5][0] = 0;face[5][1] = 3;face[5][2] = 8;face[5][3] = 0;
		face[6][0] = 7;face[6][1] = 8;face[6][2] = 3;face[6][3] = 7;
		face[7][0] = 7;face[7][1] = 4;face[7][2] = 8;face[7][3] = 7;
		face[8][0] = 5;face[8][1] = 9;face[8][2] = 1;face[8][3] = 5;
		face[9][0] = 1;face[9][1] = 9;face[9][2] = 2;face[9][3] = 1;
		face[10][0] = 6;face[10][1] = 9;face[10][2] = 5;face[10][3] = 6;
		face[11][0] = 6;face[11][1] = 2;face[11][2] = 9;face[11][3] = 6;
	}
	else{
		face[4][0] = 0;face[4][1] = 3;face[4][2] = 7;face[4][3] = 4;
		face[5][0] = 5;face[5][1] = 6;face[5][2] = 2;face[5][3] = 1;
	}
}

void set_edge_faces(){
        edge_faces[0][0] = 0;edge_faces[0][1] = 3;edge_faces[2][0] = 0;edge_faces[2][1] = 1;
        edge_faces[6][0] = 3;edge_faces[6][1] = 4;edge_faces[8][0] = 2;edge_faces[8][1] = 3;
        edge_faces[10][0] = 1;edge_faces[10][1] = 2;
        if(NEDGE == 20){
                edge_faces[1][0] = 0;edge_faces[1][1] = 9;edge_faces[3][0] = 0;edge_faces[3][1] = 5;
                edge_faces[4][0] = 1;edge_faces[4][1] = 11;edge_faces[5][0] = 1;edge_faces[5][1] = 6;
                edge_faces[7][0] = 3;edge_faces[7][1] = 8;edge_faces[9][0] = 2;edge_faces[9][1] = 7;
                edge_faces[11][0] = 2;edge_faces[11][1] = 10;edge_faces[12][0] = 4;edge_faces[12][1] = 5;
                edge_faces[13][0] = 5;edge_faces[13][1] = 6;edge_faces[14][0] = 4;edge_faces[14][1] = 7;
                edge_faces[15][0] = 6;edge_faces[15][1] = 7;edge_faces[16][0] = 8;edge_faces[16][1] = 9;
                edge_faces[17][0] = 9;edge_faces[17][1] = 11;edge_faces[18][0] = 8;edge_faces[18][1] = 10;
                edge_faces[19][0] = 10;edge_faces[19][1] = 11;
        }
        else{
                edge_faces[1][0] = 0;edge_faces[1][1] = 5;edge_faces[3][0] = 0;edge_faces[3][1] = 4;
                edge_faces[4][0] = 1;edge_faces[4][1] = 5;edge_faces[5][0] = 1;edge_faces[5][1] = 4;
                edge_faces[7][0] = 3;edge_faces[7][1] = 5;edge_faces[9][0] = 2;edge_faces[9][1] = 4;
                edge_faces[11][0] = 2;edge_faces[11][1] = 5;
        }
}

