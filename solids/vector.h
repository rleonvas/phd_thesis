#include <stdio.h>
#include <stdlib.h>
#ifndef _VECTOR_H_
#define _VECTOR_H_

void cross_product(double vec_one[3], double vec_two[3], double vec_cross[3]);

double dot_product(double vec_one[], double vec_two[], int n);

double min_value(double v_0, double v_1, double v_2);

double max_value(double v_0, double v_1, double v_2);

void merge(int *a, int n, int m);

void merge_sort(int *a, int n);

#endif
