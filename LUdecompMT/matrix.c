#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <math.h>
#include "LUdecomp.h"

#define eps 1e-50


int matrixFromFile(FILE* fin, int n, double* mat, double* mat1, double* val) {
	double curr;
	int i,
	    j,
	    count = 0;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {		
			if (fscanf(fin, "%lf", &curr) == -1) return -1;	
			else {
				mat[i * n + j] = curr;
				mat1[i * n + j] = curr;
				count++;
			}
		}
		if (fscanf(fin, "%lf", &curr) == -1) return -1;
		else {
			val[i] = curr;
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


int matrixFromFormula(int n, char* formula, double* mat, double* val) {
	int i,
	    j;
	
	if (strcmp(formula, "Identity") == 0) {
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				mat[i * n + j] = (i == j);
	}
	else {
		if (strcmp(formula, "Hilbert") == 0) {
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					mat[i * n + j] = 1 / ((double)(i + j) + 1);
		}
		else {
			if (strcmp(formula, "Lehmer") == 0) {
				for (i = 0; i < n; i++)
					for (j = 0; j < n; j++) {
						if (i >= j) mat[i * n + j] = (double)(i + 1) / (double)(j + 1);
						else mat[i * n + j] = (double)(j + 1) / (double)(i + 1);
					}
			}
			else {
				if (strcmp(formula, "NoLU") == 0) {
					for (i = 0; i < n; i++) {
						for (j = 0; j < n - i - 1; j++)
							mat[i * n + j] = 6;
						for (; j < n; j++)
							mat[i * n + j] = 5;
					}
					mat[(n - 1) * n + (n - 1)] = 4;
				}
			else return -1;
			}
		}
	}

	for (i = 0; i < n; i++) {
		val[i] = 0;
		for (j = 0; j < n; j+=2)
			val[i] += mat[i * n + j];
	}
	return 0;
}


double currentTimeSec(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec / 1000000000;
}


double currentTimeSecThread(void) {
    struct timespec t;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    return (double)t.tv_sec + (double)t.tv_nsec / 1000000000;
}


int systemSolution(int n, double* mat, double* val, double* sol, int rank, int threads, int *fl) {
	int i,
	    j,
		k,
        r1,
        r2;
    double temp;

    if (fabs(mat[0]) < eps) {
        *fl = 0;
        return -1;
    }
	for (k = 1; k < n; k++)
        mat[k] /= mat[0];
	
	for (i = 1; i < n; i++) {
        r1 = (n - i) *  rank      / threads + i;
        r2 = (n - i) * (rank + 1) / threads + i;
		
        for (k = r1; k < r2; k++) {
            temp = mat[k * n + i];
            for (j = 0; j < i; j++)
                temp -= mat[j * n + i] * mat[k * n + j];
            mat[k * n + i] = temp;
        }
        synchronize(threads);
        
		r1 = (n - i - 1) *  rank      / threads + i + 1;
		r2 = (n - i - 1) * (rank + 1) / threads + i + 1;
		
        for (k = r1; k < r2; k++) {
            temp = mat[i * n + k];
            for (j = 0; j < i; j++)
                temp -= mat[i * n + j] * mat[j * n + k];
            
            if (fabs(mat[i * n + i]) < eps) {
                *fl = 0;
                return -1;
            }
            mat[i * n + k] = temp / mat[i * n + i];
		}
		synchronize(threads);
    }

	if (rank == 0) {
		for (i = 0; i < n; i++) {
			temp = val[i];
            for (j = 0; j < i; j++)
                temp -= mat[i * n + j] * sol[j];
            if (fabs(mat[i * n + i]) < eps) {
                *fl = 0;
                return -1;
            }
			sol[i] = temp / mat[i * n + i];
		}
		
		for (i = n - 1; i >= 0; i--) {
            temp = sol[i];
			for (j = i + 1; j < n; j++)
				temp -= mat[i * n + j] * sol[j];
            sol[i] = temp;
        }
    }
    synchronize(threads);
	return 0;
}


double normFormula(int n, double* a) {
	int i;
	double res = 0;
	
	for (i = 0; i < n; i++) {
		if (i % 2 == 0)
			res += (a[i] - 1) * (a[i] - 1);
		else
			res += a[i] * a[i];
	}
	return sqrt(res);
}


double normFile(int n, double* a, double* b) {
	int i;
	double res = 0;

	for (i = 0; i < n; i++)
		res += (a[i] - b[i]) * (a[i] - b[i]);
	return sqrt(res);
}


void multiplicationMatVec(int n, double* mat, double* v, double* res) {
	int i,
	    j;
	
	for (i = 0; i < n; i++) {
		res[i] = 0;
		for (j = 0; j < n; j++)
			res[i] += mat[i * n + j] * v[j];
	}
}


void synchronize(int total_threads) {
	static pthread_mutex_t mutex       = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t  condvar_in  = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t  condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in  = 0;
	static int threads_out = 0;

	pthread_mutex_lock(&mutex);

	threads_in++;
	if (threads_in >= total_threads) {
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	}
	else
		while (threads_in < total_threads)
			pthread_cond_wait(&condvar_in,&mutex);

	threads_out++;
	if (threads_out >= total_threads) {
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	}
	else
		while (threads_out < total_threads)
			pthread_cond_wait(&condvar_out,&mutex);

	pthread_mutex_unlock(&mutex);
}
