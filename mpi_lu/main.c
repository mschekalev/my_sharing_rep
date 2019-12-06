#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h> 
#include <string.h>
#include <mpi.h>
#include "LUdecomp.h"

#define eps 1e-50


int main(int argc, char **argv)
{
    int i,
        n1 = 0,
        n = 0,
        p,
	    opt,
        rank,
        fl = 0,
        fl_res = 0,
        count = 0,
        count_res = 0,
	    f_v   = 0,
	    f_s   = 0,
	    f_f   = 0,
 	    f_n   = 0,
	    f_in  = 0,
	    f_out = 0;

    FILE *fin = stdin,
         *fout;

    double *mat = NULL,
		   *mat1 = NULL,
           *val = NULL,
           *sol,
           *vec,
           time = 0.;

	struct Params pars;
    pars.size  = 0;
    pars.scale = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

	if (argc == 1) {
        if (rank == 0) printf("No arguments\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD, 1);
	}

	pars.fin_name  = (char*)malloc(100);
	pars.fout_name = (char*)malloc(100);
	pars.formula   = (char*)malloc(100);

	while ((opt = getopt(argc, argv, "vn:f:i:o:s:")) != -1) {
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
				pars.scale = strtol(optarg, &next, 4);
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
			default: {
                if (rank == 0) printf("Correct parameters:\n -n  Matrix size, positive integer\n -l  Matrix scale, positive integer\n -i  Input file, string\n -o  Output file, string\n -f  Formula, string\n -v  Verbose mode\n");
				free(pars.fin_name);
				free(pars.fout_name);
				free(pars.formula);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Abort(MPI_COMM_WORLD, 2);
			}
		}
	}


	if (f_f) {
		if (f_n && pars.size > 0) {
			strcpy(pars.fin_name, "stdin");
			fin = stdin;
			n = pars.size;
            if (n != n) {
                if (rank == 0) printf("Too large matrix size\n");
                free(pars.fin_name);
                free(pars.fout_name);
                free(pars.formula);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Abort(MPI_COMM_WORLD, 3);
            }

            n1 = n / p;
            if (rank + 1 <= n % p)
                n1++;

            mat  = (double*)malloc(n * n1 * sizeof(double));
            //mat1 = (double*)malloc(n * n1 * sizeof(double));
			val  = (double*)malloc(n * sizeof(double));

            fl = matrixFromFormula(n, pars.formula, mat, val, rank, p);
            MPI_Allreduce(&fl, &fl_res, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            if (fl_res < 0) {
                if (rank == 0) printf("Incorect formula format. The correct ones:\n  Identity\n  Hilbert\n  Lehmer\n  NoLU\n");
				free(mat);
                //free(mat1);
				free(val);
				free(pars.fin_name);
				free(pars.fout_name);
				free(pars.formula);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Abort(MPI_COMM_WORLD, 4);
            }
		}
		else {
            if (rank == 0) printf("Matrix size must be defined and be greater than 0\n");
			free(pars.fin_name);
			free(pars.fout_name);
			free(pars.formula);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, 5);
		}
	}
	else {
		if (f_in) {
			fin = fopen(pars.fin_name, "r");
			if (!fin) {
                if (rank == 0) printf("Input file not found\n");
				free(pars.fin_name);
				free(pars.fout_name);
				free(pars.formula);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Abort(MPI_COMM_WORLD, 6);
			}
		}
		else {
            if (rank == 0) printf("No formula or file name entered\n");
			free(pars.fin_name);
			free(pars.fout_name);
			free(pars.formula);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, 7);
		}
		if (fscanf(fin, "%d", &n) == -1) {
            if (rank == 0) printf("No matrix size given\n");
			free(pars.fin_name);
			free(pars.fout_name);
			free(pars.formula);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, 8);
		}
        n1 = n / p;
        if (rank + 1 <= n % p)
            n1++;

        mat  = (double*)malloc(n * n1 * sizeof(double));
        mat1 = (double*)malloc(n * n1 * sizeof(double));
        val  = (double*)malloc(n1 * sizeof(double));

        count = matrixFromFile(fin, n, mat, mat1, val, rank, p);
        MPI_Allreduce(&count, &count_res, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (count_res != n * n + n) {
            if (rank == 0) printf("Incorrect matrix format\n");
			free(mat);
            free(mat1);
			free(val);
			free(pars.fin_name);
			free(pars.fout_name);
			free(pars.formula);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, 9);
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


    if (rank == 0 && f_v == 1) {
        printf("Processes:    %d\n", p);
		if (f_f) printf("Formula:      %s\n", pars.formula);
		else printf("Formula:      no\n");
		if (f_n) printf("Matrix size:  %d\n", pars.size);
		else printf("Matrix size:  from file\n");
		if (f_s) printf("Matrix scale: %d\n", f_s);
		else printf("Matrix scale: No\n");
       		
		printf("Input file:   %s\n", pars.fin_name);
		printf("Output file:  %s\n", pars.fout_name);
		printf("Verbose mode: 1\n\n");
	}

	
    if (f_v == 1) {
        if (rank == 0) printf("Matrix:\n");
        matrixToFile(stdout, f_s, n, mat, rank, p);
        if (rank == 0) printf("\nVector:\n");
        vectorToFile(stdout, f_s, n, val, rank, p);
        if (rank == 0) printf("\n");
	}


	time = MPI_Wtime();
    fl = LUdecomposition(n, mat, eps, rank, p);

    if (fl < 0) {
        if (rank == 0) printf("Cannot use LU decomposition on this matrix\n");
		free(mat);
		free(val);
		free(pars.fin_name);
		free(pars.fout_name);
		free(pars.formula);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD, 10);
    }

	if (f_v == 1) {
        if (rank == 0) printf("LU decomposition:\n");
        matrixToFile(stdout, f_s, n, mat, rank, p);
        if (rank == 0) printf("\n");
	}

    sol = (double*)malloc(n * sizeof(double));

    fl = systemSolution(n, mat, val, sol, eps, rank, p);
	
	time = MPI_Wtime() - time;

    if (fl < 0) {
        if (rank == 0) printf("Cannot use LU decomposition on this matrix\n");
		free(mat);
		free(val);
		free(sol);
		free(pars.fin_name);
		free(pars.fout_name);
		free(pars.formula);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD, 11);
    }

    if (rank == 0) {
        fprintf(fout, "\nSolution:\n");
		for (i = 0; i < n; i++)
			fprintf(fout, "%f ", sol[i]);
        fprintf(fout, "\n\n");
    }
		

	if (f_f && rank == 0) printf("error: %f\n", normFormula(n, sol));
	if (f_in) {
        vec = (double*)malloc(n * sizeof(double));
        multiplicationMatVec(n, mat1, sol, vec, rank, p);
        if (rank == 0) printf("error: %f\n", normFile(n, vec, val));
        free(vec);
	}

    if (rank == 0) printf("\nRun time: %f seconds\n", time);

	if (f_out) fclose(fout);
	if (f_in)  fclose(fin);
    free(mat);
    free(mat1);
    free(val);
    free(sol);
    free(pars.fin_name);
    free(pars.fout_name);
    free(pars.formula);
    MPI_Finalize();
	return 0;
}
