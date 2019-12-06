#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <fenv.h>
#include "QRdecomp.h"


int matrixFromFile(FILE* fin, int n, double* mat, double* mat1) {
	double curr;
	int i,
	    j,
	    count = 0;
	
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {		
			if (fscanf(fin, "%lf", &curr) == -1) return -1;	
			else {
				mat[i * n + j] = curr;
				mat1[i * n + j] = curr;
				count++;
			}
		}
	return count;
}


void matrixToFile(FILE* fout, int s, int n, double* mat) {
	int i,
	    j,
		lim = s;

	if (s >= n - 1) lim = n;
	for (i = 0; i < lim; i++) {
		for (j = 0; j < lim; j++)
			fprintf(fout, "%f ", mat[i * n + j]);
		if (s < n - 1) fprintf(fout, ". . . %f", mat[i * n + n - 1]);
		fprintf(fout, "\n");
	}
	if (s < n - 1) {
		fprintf(fout, ". . .\n. . .\n. . .\n");
		for (j = 0; j < lim; j++)
			fprintf(fout, "%f ", mat[(n - 1) * n + j]);
		fprintf(fout, ". . . %f\n", mat[(n - 1) * n + n - 1]);
	}
}


int matrixFromFormula(int n, char* formula, double* mat, double* mat1) {
	int i,
	    j;
	
	if (strcmp(formula, "Identity") == 0)
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				mat[i * n + j]  = (i == j);
				mat1[i * n + j] = (i == j);
			}
	else {
		if (strcmp(formula, "Hilbert") == 0)
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++) {
					mat[i * n + j]  = 1 / ((double)(i + j) + 1);
					mat1[i * n + j] = 1 / ((double)(i + j) + 1);
				}
		else {
			if (strcmp(formula, "Lehmer") == 0)
				for (i = 0; i < n; i++)
					for (j = 0; j < n; j++) {
						if (i >= j) {
							mat[i * n + j]  = (double)(i + 1) / (double)(j + 1);
							mat1[i * n + j] = (double)(i + 1) / (double)(j + 1);
						}	
						else {
							mat[i * n + j]  = (double)(j + 1) / (double)(i + 1);
							mat1[i * n + j] = (double)(j + 1) / (double)(i + 1);
						}
					}
			else {
				if (strcmp(formula, "NoLU") == 0) {
					for (i = 0; i < n; i++) {
						for (j = 0; j < n - i - 1; j++) {
							mat[i * n + j]  = 6;
							mat1[i * n + j] = 6;
						}
						for (; j < n; j++) {
							mat[i * n + j]  = 5;
							mat1[i * n + j] = 5;
						}
					}
					mat[(n - 1) * n + (n - 1)]  = 4;
					mat1[(n - 1) * n + (n - 1)] = 4;
				}
				else {
					if (strcmp(formula, "Sym") == 0)
						for (i = 0; i < n; i++)
							for (j = 0; j < n; j++) {
								mat[i * n + j]  = fabs(i - j);
								mat1[i * n + j] = fabs(i - j);
							}
                    else {
                        if (strcmp(formula, "2.10") == 0) {
                            for (i = 0; i < n - 1; i++)
                                for (j = 0; j < n - 1; j++) {
                                    if (i == j) {
                                        mat[i * n + j]  = 1;
                                        mat1[i * n + j] = 1;
                                    }
                                    else {
                                        mat[i * n + j]  = 0;
                                        mat1[i * n + j] = 0;
                                    }
                                }
                            for (i = 0; i < n; i++) {
                                mat[i * n + n - 1]    = i + 1;
                                mat[(n - 1) * n + i]  = i + 1;
                                mat1[i * n + n - 1]   = i + 1;
                                mat1[(n - 1) * n + i] = i + 1;
                            }
                        }
						else {
							if (strcmp(formula, "2.15") == 0) {
								for (i = 0; i < n; i++) {
									for (j = i + 1; j < n; j++) {
										mat[i * n + j]  = 2;
										mat1[i * n + j] = 2;
									}
									mat[i * n + i]  = i + 2;
									mat1[i * n + i] = i + 2;
									mat[i * n]  = -(n - i);
									mat1[i * n] = -(n - i);
									for (j = 1; j < i; j++) {
										mat[i * n + j]  = 0;
										mat1[i * n + j] = 0;
									}
								}
								mat[0] -= n;
							}
							else return -1;
						}
                    }
				}
			}
		}
	}
	return 0;
}


/*double currentTimeSec(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec / 1000000000;
}*/


void toAlmostTriangular(int n, double* mat, double eps) {
	int i,
		j,
		k;
	double cos,
		   sin,
		   temp1,
           temp2,
		   sq;
	
	for (i = 1; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
            temp1 = mat[i * n + i - 1];
            temp2 = mat[j * n + i - 1];
			sq = sqrt(temp1 * temp1 + temp2 * temp2);
			if (sq > 0) {
				cos =  temp1 / sq;
				sin = -temp2 / sq;
				
                for (k = i - 1; k < n; k++) {
					temp1 = mat[i * n + k];
                    temp2 = mat[j * n + k];
					mat[i * n + k] = temp1 * cos - temp2 * sin;
					mat[j * n + k] = temp1 * sin + temp2 * cos;
				}
				
                for (k = 0; k < n; k++) {
					temp1 = mat[k * n + i];
                    temp2 = mat[k * n + j];
					mat[k * n + i] = temp1 * cos - temp2 * sin;
					mat[k * n + j] = temp1 * sin + temp2 * cos;
				}
			}
		}
	}
}


void qrForAlmostTriangular(int n, int acs, double* mat, double* vcos, double* vsin, double eps) {
	int i,
		j;
    double temp1,
           temp2,
		   sq;
	
	for (i = 0; i <= acs - 1; i++) {
        temp1 = mat[i * n + i];
        temp2 = mat[(i + 1) * n + i];
		sq = sqrt(temp1 * temp1 + temp2 * temp2);
		if (sq > 0) {
			vcos[i] =  temp1 / sq;
			vsin[i] = -temp2 / sq;
			
            for (j = i; j <= acs; j++) {
				temp1 = mat[i * n + j];
                temp2 = mat[(i + 1) * n + j];
				mat[i * n + j]       = temp1 * vcos[i] - temp2 * vsin[i];
				mat[(i + 1) * n + j] = temp1 * vsin[i] + temp2 * vcos[i];
			}
		}
		else {
			vcos[i] = 1;
			vsin[i] = 0;
		}
	}
}


void rqMultiplication(int n, int acs, double* mat, double* vcos, double* vsin) {
	int i,
		j;
	double temp1,
           temp2;
	
	for (i = 0; i <= acs - 1; i++)
        for (j = 0; j < i + 2; j++) {
			temp1 = mat[j * n + i];
            temp2 = mat[j * n + i + 1];
			mat[j * n + i]     = temp1 * vcos[i] - temp2 * vsin[i];
			mat[j * n + i + 1] = temp1 * vsin[i] + temp2 * vcos[i];
		}
}


void matrixShift(int n, int acs, double* mat, double sh) {
	int i;
	
	for (i = 0; i <= acs; i++)
		mat[i * n + i] -= sh;
}


double qrIteration(int n, int acs, double* mat, double* vcos, double* vsin, double eps) {
	int i,
		j;
	double val,
		   sh;
	
	while (fabs(mat[acs * n + acs - 1]) > eps) {
		sh = mat[acs * n + acs];
		
		matrixShift(n, acs, mat, sh);
		qrForAlmostTriangular(n, acs, mat, vcos, vsin, eps);
		rqMultiplication(n, acs, mat, vcos, vsin);
		matrixShift(n, acs, mat, -sh);
	}
	val = mat[acs * n + acs];
	return val;
}


void qrAlgorithm(int n, double* mat, double* vcos, double* vsin, double* sol, double eps) {
	int i;
	double a,
		   b,
		   c,
		   d;
	
	for (i = n - 1; i >= 2; i--)
		sol[i] = qrIteration(n, i, mat, vcos, vsin, eps);
	
	a = mat[0];
	b = mat[1];
	c = mat[n];
	d = mat[n + 1];
	sol[0] = ((a + d) + sqrt((a + d) * (a + d) - 4 * (a * d - c * b))) / 2;
	sol[1] = ((a + d) - sqrt((a + d) * (a + d) - 4 * (a * d - c * b))) / 2;
}


double normTrace(int n, double* mat, double* sol) {
    int i;
    double res = 0;

    for (i = 0; i < n; i++)
        res += mat[i * n + i] - sol[i];
    return res;
}


double normMatrix(int n, double* mat, double* sol) {
    int i,
        j;
    double res1 = 0,
           res2 = 0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            res1 += (mat[i * n + j] * mat[i * n + j]);
        res2 += (sol[i] * sol[i]);
    }
    return sqrt(res1) - sqrt(res2);
}


void vectorToFile(FILE* fout, int s, int n, double* vec) {
	int i,
		lim = s;
	
	if (s >= n - 1) lim = n;
	for (i = 0; i < lim; i++)
		fprintf(fout, "%f ", vec[i]);
    if (s < n - 1) fprintf(fout, " . . . %f", vec[n - 1]);
	fprintf(fout, "\n");
}
