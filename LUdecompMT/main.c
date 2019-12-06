#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h> 
#include <string.h>
#include "LUdecomp.h"


double thread_time = 0;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;


void* Solution(void *pa) {
    Args* arg = (Args*)pa;
    double time = 0;

    time = currentTimeSecThread();
    systemSolution(arg->n, arg->mat, arg->val, arg->sol, arg->rank, arg->threads, &(arg->fl));
    time =  currentTimeSecThread()- time;

    arg->time = time;

    pthread_mutex_lock(&mutex);
    thread_time += time;
    pthread_mutex_unlock(&mutex);
    return NULL;
}


int main(int argc, char** argv)
{
	
	int i,
	    n,
	    opt,
        f_t   = 0,
	    f_v   = 0,
	    f_s   = 0,
	    f_f   = 0,
 	    f_n   = 0,
	    f_in  = 0,
	    f_out = 0;

	FILE* fin,
	    * fout;

	double* mat,
	      * mat1,
	      * val,
	      * sol,
          * vect,
            time = 0;

    struct Params pars;
    Args*      args;
	pthread_t* threads;

	if (argc == 1) {
		printf("No arguments\n");
		return -1;
	}
	pars.fin_name  = (char*)malloc(100);
	pars.fout_name = (char*)malloc(100);
	pars.formula   = (char*)malloc(100);

	while ((opt = getopt(argc, argv, "vn:f:i:o:s:t:")) != -1) {
		switch (opt) {
			case 'v': {
		        pars.verbose = 1;
				f_v = 1;
		        break;
		    	}
			case 'n': {
				char* next;
				pars.size = strtol(optarg, &next, 10);
				f_n = 1;
		        break;
			}
			case 's': {
				char* next;
				pars.scale = strtol(optarg, &next, 10);
				f_s = pars.scale;
				break;
			}
            case 't': {
				char* next;
				pars.threads = strtol(optarg, &next, 10);
				f_t = pars.threads;
				break;
			}
			case 'i': {
				strcpy(pars.fin_name, optarg);
				f_in = 1;
				break;
			}
			case 'o': {
				strcpy(pars.fout_name, optarg);
				f_out = 1;
				break;
			}
			case 'f': {
				strcpy(pars.formula, optarg);
				f_f = 1;
				break;
			}
			default: {
				printf("Correct parameters:\n -n  Matrix size, positive integer\n -l  Matrix scale, positive integer\n -t  Threads number, positive integer\n -i  Input file, string\n -o  Output file, string\n -f  Formula, string\n -v  Verbose mode\n");
                free(pars.fin_name);
                free(pars.fout_name);
                free(pars.formula);
				return -2;
			}
		}
	}


	if (f_f) {
		if (f_n && pars.size > 0) {
			strcpy(pars.fin_name, "stdin");
			fin = stdin;
			n = pars.size;
			mat  = (double*)malloc(n * n * sizeof(double));
			mat1 = (double*)malloc(n * n * sizeof(double));
			val  = (double*)malloc(n * sizeof(double));
			if (matrixFromFormula(n, pars.formula, mat, val) != 0) {
				printf("Incorect formula format. The correct ones:\n  Identity\n  Hilbert\n  Lehmer\n  NoLU\n");
                free(mat);
                free(val);
                free(pars.fin_name);
                free(pars.fout_name);
                free(pars.formula);
				return -3;
			}
		}
		else {
			printf("Matrix size must be defined and be greater than 0\n");
            free(pars.fin_name);
            free(pars.fout_name);
            free(pars.formula);
			return -4;
		}
	}
	else {
		if (f_in) {
			fin = fopen(pars.fin_name, "r");
			if (!fin) {
				printf("Input file not found\n");
                free(pars.fin_name);
                free(pars.fout_name);
                free(pars.formula);
				return -5;
			}
		}
		else {
			printf("No formula or file name entered\n");
            free(pars.fin_name);
            free(pars.fout_name);
            free(pars.formula);
			return -6;
		}
		if (fscanf(fin, "%d", &n) == -1) {
			printf("No matrix size given\n");
            free(pars.fin_name);
            free(pars.fout_name);
            free(pars.formula);
			return -7;
		}
        mat  = (double*)malloc(n * n * sizeof(double));
		mat1 = (double*)malloc(n * n * sizeof(double));
		val  = (double*)malloc(n * sizeof(double));
		if (matrixFromFile(fin, n, mat, mat1, val) != n * n + n) {
			printf("Incorrect matrix format\n");
            free(mat);
            free(mat1);
            free(val);
            free(pars.fin_name);
            free(pars.fout_name);
            free(pars.formula);
			return -8;
		}
	}
	if (!f_out) {
		strcpy(pars.fout_name, "stdout");
		fout = stdout;
	}
	else
		fout = fopen(pars.fout_name, "w");
	if (!f_s || pars.scale > n || pars.scale <= 0) {
		pars.scale = n;
		f_s = n;
	}
    if (f_t <= 0) f_t = 4;
	
	if (f_v == 1) {
		if (f_f) printf("Formula:      %s\n", pars.formula);
		else printf("Formula:      no\n");
		if (f_n) printf("Matrix size:  %d\n", pars.size);
		else printf("Matrix size:  from file\n");
		if (f_s) printf("Matrix scale: %d\n", f_s);
		else printf("Matrix scale: No\n");
        printf("Threads:      %d\n", f_t);
		printf("Input file:   %s\n", pars.fin_name);
		printf("Output file:  %s\n", pars.fout_name);
		printf("Verbose mode: 1\n\n");
	}
	
	threads = (pthread_t*)malloc(f_t * sizeof(pthread_t));
    args    = (Args*)malloc(f_t * sizeof(Args));
	sol     = (double*)malloc(n * sizeof(double));
	
	if (f_v == 1) {
		printf("Matrix:\n");
		matrixToFile(stdout, f_s, n, mat);
		printf("\n");
	}
	
	for (i = 0; i < f_t; i++) {
        args[i].n       = n;
        args[i].fl      = 1;
        args[i].mat     = mat;
        args[i].val     = val;
        args[i].sol     = sol;
        args[i].rank    = i;
		args[i].threads = f_t;
        args[i].time    = 0;
	}

    time = currentTimeSec();
    for (i = 0; i < f_t; i++)
        if (pthread_create(threads + i, 0, Solution, args + i)) {
            printf("Cannot create thread  %d!\n", i);
            free(mat);
            free(mat1);
            free(val);
            free(sol);
            free(args);
            free(threads);
            free(pars.fin_name);
            free(pars.fout_name);
            free(pars.formula);
            return -9;
		}

    for (i = 0; i < f_t; i++) {
		if (pthread_join(threads[i], 0)) {
            printf("Cannot wait thread  %d!\n", i);
            free(mat);
            free(mat1);
            free(val);
            free(sol);
            free(args);
            free(threads);
            free(pars.fin_name);
            free(pars.fout_name);
            free(pars.formula);
            return -10;
        }
        if (args[i].fl == 0) {
            printf("LU decomposition cannot be done\n");
            free(mat);
            free(mat1);
            free(val);
            free(sol);
            free(args);
            free(threads);
            free(pars.fin_name);
            free(pars.fout_name);
            free(pars.formula);
            return -11;
        }
    }
    time = currentTimeSec() - time;
	
	fprintf(fout, "\nSolution:\n");
	for (i = 0; i < n; i++)
		fprintf(fout, "%f ", sol[i]);
    fprintf(fout, "\n");
    if (f_f) printf("\nError: %f\n", normFormula(n, sol));
	if (f_in) {
		vect = (double*)malloc(n * sizeof(double));
		multiplicationMatVec(n, mat1, sol, vect);
        printf("\nError: %f\n", normFile(n, vect, val));
		free(vect);
	}

    printf("\nRun time:          %f seconds", time);
    printf("\nThread run time:   %f seconds\n\n", thread_time);
    for (i = 0; i < f_t; i++)
        printf("Thread %d run time: %f seconds\n", i, args[i].time);
    //printf("\nPer thread time: %f seconds\n", thread_time / f_t);
	if (f_out) fclose(fout);
    if (f_in)  fclose(fin);
    free(mat);
    free(mat1);
    free(val);
    free(sol);
    free(args);
    free(threads);
    free(pars.fin_name);
    free(pars.fout_name);
    free(pars.formula);
	return 0;
}
