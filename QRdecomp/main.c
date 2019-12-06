#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h> 
#include <string.h>
#include <fenv.h>
#include "QRdecomp.h"


int main(int argc, char** argv)
{
    //int i,
    int n = 5,
	    opt,
	    f_v   = 0,
	    f_s   = 0,
	    f_f   = 0,
 	    f_n   = 0,
	    f_in  = 0,
        f_e   = 0,
	    f_out = 0;

	FILE* fin,
	    * fout;

	double* mat,
          * mat1,
          * sol,
		  * vcos,
          * vsin,
            eps = -1;
            //time;

	struct Params pars;

	if (argc == 1) {
		printf("No arguments\n");
		return -7;
	}
	pars.fin_name  = (char*)malloc(100);
	pars.fout_name = (char*)malloc(100);
	pars.formula   = (char*)malloc(100);

	while ((opt = getopt(argc, argv, "vn:f:i:o:s:e:")) != -1) {
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
            case 'e': {
                char* next;
				pars.eps = strtod(optarg, &next);
                eps = pars.eps;
                f_e = 1;
		        break;
            }
			default: {
				printf("Correct parameters:\n -n  Matrix size, positive integer\n -l  Matrix scale, positive integer\n -i  Input file, string\n -o  Output file, string\n -f  Formula, string\n -v  Verbose mode\n -e  Epsilon\n");
				free(pars.fin_name);
				free(pars.fout_name);
				free(pars.formula);
				return -6;
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
			if (matrixFromFormula(n, pars.formula, mat, mat1)) {
                printf("Incorect formula format. The correct ones:\n  Identity\n  Hilbert\n  Lehmer\n  NoLU\n  Sym\n  2.10\n  2.15\n");
				free(mat);
                free(mat1);
				free(pars.fin_name);
				free(pars.fout_name);
				free(pars.formula);
				return -5;
			}
		}
		else {
			printf("Matrix size must be defined and be greater than 0\n");
			free(pars.fin_name);
			free(pars.fout_name);
			free(pars.formula);
			return -1;
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
				return -2;
			}
		}
		else {
			printf("No formula or file name entered\n");
			free(pars.fin_name);
			free(pars.fout_name);
			free(pars.formula);
			return -8;
		}
		if (fscanf(fin, "%d", &n) == -1) {
			printf("No matrix size given\n");
			free(pars.fin_name);
			free(pars.fout_name);
			free(pars.formula);
			return -3;
		}
        mat  = (double*)malloc(n * n * sizeof(double));
        mat1 = (double*)malloc(n * n * sizeof(double));
		if (matrixFromFile(fin, n, mat, mat1) != n * n) {
			printf("Incorrect matrix format\n");
			free(mat);
            free(mat1);
			free(pars.fin_name);
			free(pars.fout_name);
			free(pars.formula);
			return -4;
		}
	}
	if (!f_out) {
		strcpy(pars.fout_name, "stdout");
		fout = stdout;
	}
	else
		fout = fopen(pars.fout_name, "w");
    if (!f_s || pars.scale > n) {
        pars.scale = n;
        f_s = n;
	}
    if (!f_e) {
        pars.eps = 1e-6;
        eps = 1e-6;
    }

	if (f_v == 1) {
		if (f_f) printf("Formula:      %s\n", pars.formula);
		else printf("Formula:      no\n");
        if (f_n) printf("Matrix size:  %d\n", n);
		else printf("Matrix size:  from file\n");
		if (f_s) printf("Matrix scale: %d\n", f_s);
		else printf("Matrix scale: no\n");
        printf("Precision:    %.12f\n", eps);
        printf("Input file:   %s\n", pars.fin_name);
        printf("Output file:  %s\n", pars.fout_name);
		printf("Verbose mode: 1\n\n");
	}

	
	if (f_v == 1) {
		fprintf(stdout, "Matrix:\n");
		matrixToFile(stdout, f_s, n, mat1);
		fprintf(stdout, "\n");
	}
    //time = currentTimeSec();
	
    toAlmostTriangular(n, mat, eps);
	if (f_v == 1) {
		fprintf(stdout, "Almost triangular matrix:\n");
		matrixToFile(stdout, f_s, n, mat);
		fprintf(stdout, "\n");
	}
	
	vcos = (double*)malloc((n - 1) * sizeof(double));
	vsin = (double*)malloc((n - 1) * sizeof(double));
    sol  = (double*)malloc(n * sizeof(double));
	
    qrAlgorithm(n, mat, vcos, vsin, sol, eps);
	
	if (f_v == 1) {
		fprintf(stdout, "Matrix after QR-algorithm:\n");
		matrixToFile(stdout, f_s, n, mat);
        fprintf(stdout, "\n\n");
	}

    //time = currentTimeSec() - time;

    fprintf(stdout, "Eigen values:\n");
    vectorToFile(fout, n, n, sol);
	
    fprintf(stdout, "\nTrace difference:    %f",   normTrace(n, mat1, sol));
    fprintf(stdout, "\nVector difference:   %f\n", normMatrix(n, mat1, sol));
	
    //fprintf(stdout, "\nRun time: %f seconds\n", time);
	if (f_out) fclose(fout);
	if (f_in)  fclose(fin);
	free(mat);
    free(mat1);
    free(sol);
    free(vcos);
    free(vsin);
	free(pars.fin_name);
	free(pars.fout_name);
	free(pars.formula);
	return 0;
}
