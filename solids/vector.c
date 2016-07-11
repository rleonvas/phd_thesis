#include "vector.h"

void cross_product(double vec_one[3], double vec_two[3], double vec_cross[3]){
/*Input:vec_one: vector of dimension 3
		vec_two: vector of dimension 3
		vec_cross: vector of dimension 3, where it stores the cross product result.
Output:	None.
Description: Computing cross product between two vectors of dimension 3.*/

	vec_cross[0] = vec_one[1]*vec_two[2] - vec_two[1]*vec_one[2];
	vec_cross[1] = vec_one[2]*vec_two[0] - vec_two[2]*vec_one[0];
	vec_cross[2] = vec_one[0]*vec_two[1] - vec_two[0]*vec_one[1];
}

double dot_product(double vec_one[], double vec_two[], int n){
/*Input:vec_one: vector of dimension n
		vec_two: vector of dimension n
		n: dimension of vectors.
Output:	scalar.
Description: Computing dot product between two vectors of dimension n.*/
	double dot_p = 0.0;
	int i;

	for(i = 0; i < n; i++)
		dot_p += vec_one[i]*vec_two[i];
	return dot_p;
}

double min_value(double v_0, double v_1, double v_2){
	if (v_0 < v_1 && v_0 < v_2)
			return v_0;
	else if (v_1 < v_0 && v_1 < v_2)
			return v_1;
	else
		return v_2;
}

double max_value(double v_0, double v_1, double v_2){
	if (v_0 > v_1 && v_0 > v_2)
			return v_0;
	else if (v_1 > v_0 && v_1 > v_2)
			return v_1;
	else
		return v_2;
}

void merge(int *a, int n, int m){
  int i, j, k;
  int *x = malloc(n*sizeof(int));

  for (i = 0, j = m, k = 0; k < n; k++) {
    x[k] = j == n      ? a[i++]
         : i == m      ? a[j++]
         : a[j] < a[i] ? a[j++]
         :               a[i++];
  }
  for (i = 0; i < n; i++) {
    a[i] = x[i];
  }
  free(x);
}

void merge_sort(int *a, int n){
  if (n < 2)
      return;
  int m = n / 2;
  merge_sort(a, m);
  merge_sort(a + m, n - m);
  merge(a, n, m);
}
