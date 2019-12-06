#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include "LUdecomp.h"


int matrixFromFile(FILE* fin, int n, double* mat, double* mat1, double* val, int rank, int p) {
	int i,
        j,
        count = 0;
    double curr,
          *tmp_mff = (double*)malloc((n + 1) * sizeof(double));
	MPI_Status status;

    for (i = 0; i < n; i++) {
		if (rank == 0) {
            for (j = 0; j < n + 1; j++) {
                if (fscanf(fin, "%lf", &curr) == -1) {
                    free(tmp_mff);
                    return count;
                }

                tmp_mff[j] = curr;
                count++;
			}
            if (i % p != 0) {
                MPI_Send(tmp_mff, n, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD);
				MPI_Send(tmp_mff, n, MPI_DOUBLE, i % p, 1, MPI_COMM_WORLD);
                MPI_Send(tmp_mff + n, 1, MPI_DOUBLE, i % p, 2, MPI_COMM_WORLD);
            }
            else {
                for(j = 0; j < n; j++) {
                    mat [(i / p) * n + j] = tmp_mff[j];
					mat1[(i / p) * n + j] = tmp_mff[j];
				}
                val[i] = tmp_mff[n];
            }
		}
		else
            if (rank == i % p) {
                MPI_Recv(mat  + (i / p) * n, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(mat1 + (i / p) * n, n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(val  + (i / p), 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
            }
	}
    free(tmp_mff);
    return count;
}


void matrixToFile(FILE* fout, int s, int n, double* mat, int rank, int p) {
	int i,
        j,
        lim = s;
    double *tmp_mtf = (double*)malloc(n * n * sizeof(double));
	MPI_Status status;
	
    if (s >= n - 1) lim = n;
	
    if (rank == 0) {
        for (i = 0; i < n; i++) {
			if (i % p == 0)
                for (j = 0; j < n; j++)
                    tmp_mtf[i * n + j] = mat[(i / p) * n + j];

            else
                MPI_Recv(tmp_mtf + i * n, n, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &status);
		}
        for (i = 0; i < lim; i++) {
            for (j = 0; j < lim; j++)
                fprintf (fout, "%f ", tmp_mtf[i * n + j]);
            if (s < n - 1) fprintf(fout, ". . . %f", tmp_mtf[i * n + n - 1]);
            fprintf(fout, "\n");
		}
        if (s < n - 1) {
            fprintf(fout, ". . .\n. . .\n. . .\n");
            for (j = 0; j < lim; j++)
                fprintf(fout, "%f ", tmp_mtf[(n - 1) * n + j]);
            fprintf(fout, ". . . %f\n", tmp_mtf[(n - 1) * n + n - 1]);
		}
	}
	else
        for (i = 0; i < n; i++)
            if (rank == i % p)
                MPI_Send(mat + (i / p) * n, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    free(tmp_mtf);
}


int matrixFromFormula(int n, char* formula, double* mat, double* val, int rank, int p) {
	int i,
        j;
		
    if (strcmp(formula, "Identity") == 0) {
        for (i = 0; i < n; i++)
            if (rank == i % p)
				for (j = 0; j < n; j++)
					mat[(i / p) * n + j] = (i == j);
    } else {
        if (strcmp(formula, "Hilbert") == 0) {
            for (i = 0; i < n; i++)
                if (rank == i % p)
					for (j = 0; j < n; j++)
						mat[(i / p) * n + j] = 1 / ((double)(i + j) + 1);
        } else {
            if (strcmp(formula, "Lehmer") == 0) {
                for (i = 0; i < n; i++)
                    if (rank == i % p)
						for (j = 0; j < n; j++) {
							if (i >= j) 
								mat[(i / p) * n + j] = (double)(i + 1) / (double)(j + 1);
							else 
								mat[(i / p) * n + j] = (double)(j + 1) / (double)(i + 1);
						}
            } else {
				if (strcmp(formula, "NoLU") == 0) {
                    for (i = 0; i < n; i++)
                        if (rank == i % p) {
							for (j = 0; j < n - i - 1; j++)
                                mat[(i / p) * n + j] = 6.;
							for (; j < n; j++)
                                mat[(i / p) * n + j] = 5.;
						}
                    if (rank == (n - 1) % p)
                        mat[((n - 1) / p) * n + (n - 1)] = 4.;
                } else {
                    if (strcmp(formula, "Zero") == 0) {
                        for (i = 0; i < n; i++)
                            if (rank == i % p)
                                for (j = 0; j < n; j++) {
                                    if (i == n - 1)
                                        mat[(i / p) * n + j] = 0.;
                                    else
                                        mat[(i / p) * n + j] = (i == j);
                                }
                    }
                    else
                        return -1;
                 }
            }
		}
    }

    for (i = 0; i < n; i++) {
        if (rank == i % p) {
            val[i / p] = 0.;
			for (j = 0; j < n; j+=2)
                val[i / p] += mat[(i / p) * n + j];
        }
    }
    return 0;
}


void vectorToFile(FILE *fout, int s, int n, double* val, int rank, int p) {
    int i,
        lim = s;
    double *tmp_vtf = (double*)malloc(n * sizeof(double));
    MPI_Status status;

    if (s >= n - 1) lim = n;

    for (i = 0; i < n; i++) {
        if (i % p == 0)
            tmp_vtf[i] = val[i / p];
        else {
            if (rank == 0)
                MPI_Recv(tmp_vtf + i, 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &status);

            if (rank == i % p)
                MPI_Send(val + i / p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        for (i = 0; i < lim; i++)
            fprintf(fout, "%f ", tmp_vtf[i]);
        if (s < n - 1) fprintf(fout, ". . . %f", tmp_vtf[n - 1]);
        fprintf(fout, "\n");
    }
    free(tmp_vtf);
}


int LUdecomposition(int n, double *mat, double eps, int rank, int p) {
    int i,
        k,
        j,
        l,
        n1 = (n - 1) / p + 1,
		fl = 0,
       *recvcs,
       *displs;
    double tmp  = 0.,
		   time1  = 0., time2  = 0., time3  = 0., time4  = 0.,
		   time1r = 0., time2r = 0., time3r = 0., time4r = 0.,
		  *buf,
          *buf2,
          *buf1,
          *buf3;
			
    MPI_Status status;
	

	if (rank == 0) {
        tmp = mat[0];
        if (fabs(tmp) < eps)
		    fl = -1;
	}
    MPI_Bcast(&fl, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (fl == -1) return -1;


    buf  = (double*)malloc((n - 1) * sizeof(double));
    buf1 = (double*)malloc(n1 * sizeof(double));
    buf2 = (double*)malloc(n1 * sizeof(double));
    buf3 = (double*)malloc((n - 1) * sizeof(double));
    recvcs = (int*)malloc(p * sizeof(int));
    displs = (int*)malloc(p * sizeof(int));


    if (rank == 0)
        for (k = 1; k < n; k++)
            mat[k] /= tmp;

    for (i = 1; i < n; i++) {

        time1 = MPI_Wtime();

        n1 = (i - 1) / p;
        if (rank <= (i - 1) % p)
            n1++;

        for (j = 0; j < n1; j++)
            buf1[j] = mat[j * n + i];
		
        l = 0;
        for (k = 0; k < p; k++) {
            recvcs[k] = (i - 1) / p;
            if (k <= (i - 1) % p)
                recvcs[k]++;

            displs[k] = l;
            l += recvcs[k];
        }

        MPI_Allgatherv(buf1, n1, MPI_DOUBLE, buf, recvcs, displs, MPI_DOUBLE, MPI_COMM_WORLD);


        time1 = MPI_Wtime() - time1;
        time2 = MPI_Wtime();


        for (k = i; k < n; k++) {
            
            if (rank == k % p) {
                tmp = mat[(k / p) * n + i];
                for (l = 0; l < p; l++)
                    for (j = 0; j < recvcs[l]; j++)
                        tmp -= mat[(k / p) * n + l + j * p] * buf[displs[l] + j];
                mat[(k / p) * n + i] = tmp;
            }
        }

        time2 = MPI_Wtime() - time2;

        MPI_Barrier(MPI_COMM_WORLD);

        time3 = MPI_Wtime();

        if (rank == i % p)
			for (k = 0; k < p; k++) {
				
			    if (k != i % p) {
                    l = (i - 1) / p;
				    if (k <= (i - 1) % p)
					    l++;

				    for (j = 0; j < l; j++)
                        buf2[j] = mat[(i / p) * n + k + j * p];

                    MPI_Send(buf2, l, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
                }
                else
                    for (j = 0; j < n1; j++)
                        buf1[j] = mat[(i / p) * n + k + j * p];
			}
		else
			MPI_Recv(buf1, n1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &status);
		
		
        time3 = MPI_Wtime() - time3;
        time4 = MPI_Wtime();
		

        for (k = i + 1; k < n; k++) {

            buf3[k] = 0.;
            for (j = 0; j < n1; j++)
                buf3[k] -= buf1[j] * mat[j * n + k];

            if (rank == i % p)
                buf3[k] += mat[(i / p) * n + k];
        }

        MPI_Reduce(buf3 + i + 1, mat + (i / p) * n + i + 1, n - i - 1, MPI_DOUBLE, MPI_SUM, i % p, MPI_COMM_WORLD);

        if (rank == i % p) {
            tmp = mat[(i / p) * n + i];

            if (fabs(tmp) < eps)
                fl = -1;
            else
                for (k = i + 1; k < n; k++)
                    mat[(i / p) * n + k] /= tmp;
        }
        MPI_Bcast(&fl, 1, MPI_INT, i % p, MPI_COMM_WORLD);
        if (fl == -1) {
            free(recvcs);
            free(displs);
            free(buf);
            free(buf1);
            free(buf2);
            //free(buf3);
            return -1;
        }

        time4 = MPI_Wtime() - time4;

        time1r += time1;
        time2r += time2;
        time3r += time3;
        time4r += time4;

        MPI_Barrier(MPI_COMM_WORLD);
		
	}

    if (rank == 0) {
        printf("BLOCK 1:   %f\n", time1r);
        printf("BLOCK 2:   %f\n", time2r);
        printf("BLOCK 3:   %f\n", time3r);
        printf("BLOCK 4:   %f\n\n", time4r);
    }

    free(recvcs);
    free(displs);
    free(buf);
    free(buf1);
    free(buf2);
    //free(buf3);
    return 0;
}


int systemSolution(int n, double* mat, double* val, double *sol, double eps, int rank, int p) {
    int i,
        j,
        fl = 0;
    double tmp = 0.;

    for (i = 0; i < n; i++) {

		if (rank == i % p) {
			
            tmp = val[i / p];
			for (j = 0; j < i; j++)
                tmp -= mat[(i / p) * n + j] * sol[j];
			
			if (fabs(mat[(i / p) * n + i]) < eps) 
				fl = -1;
			else
                sol[i] = tmp / mat[(i / p) * n + i];
		}
		
        MPI_Bcast(&fl, 1, MPI_INT, i % p, MPI_COMM_WORLD);
		if (fl == -1) return -1;
		
        MPI_Bcast(sol + i, 1, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
	}
	
	
	for (i = n - 1; i >= 0; i--) {
		
		if (rank == i % p) {
			
            tmp = sol[i];
			for (j = i + 1; j < n; j++)
                tmp -= mat[(i / p) * n + j] * sol[j];
            sol[i] = tmp;
		}
		
        MPI_Bcast(sol + i, 1, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
	}

    return 0;
}


double normFormula(int n, double* a) {
	int i;
	double res = 0.,
		   tmp = 0.;
	
	for (i = 0; i < n; i++) {

		tmp = a[i];
		if (i % 2 == 0)
			res += (tmp - 1) * (tmp - 1);
		else
			res += tmp * tmp;
	}
	return sqrt(res);
}


double normFile(int n, double* a, double* b) {
	int i;
	double res  = 0.,
		   tmpa = 0.,
		   tmpb = 0.;

	for (i = 0; i < n; i++) {
		tmpa = a[i];
		tmpb = b[i];
		res += (tmpa - tmpb) * (tmpa - tmpb);
	}
	return sqrt(res);
}


void multiplicationMatVec(int n, double* mat, double* v, double* res, int rank, int p) {
	int i,
	    j;
	
	for (i = 0; i < n; i++) {
		res[i] = 0.;
		if (rank == i % p)
			for (j = 0; j < n; j++)
				res[i / p] += mat[(i / p) * n + j] * v[j];
	}
}
