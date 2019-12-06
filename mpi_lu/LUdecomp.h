struct Params {
 	int   size;
	int   scale;
	char  verbose;
	char* fin_name;
	char* fout_name;
	char* formula;
};


double normFormula(int n, double* a);
double normFile(int n, double* a, double* b);
void   multiplicationMatVec(int n, double* mat, double* v, double* res, int rank, int p);
void   matrixToFile(FILE* fout, int s, int n, double* mat, int rank, int p);
int    systemSolution(int n, double* mat, double* val, double *buf, double eps, int rank, int p);
int    LUdecomposition(int n, double *mat, double eps, int rank, int p);
int    matrixFromFile(FILE* fin, int n, double* mat, double* mat1, double* val, int rank, int p);
int    matrixFromFormula(int n, char* formula, double* mat, double* val, int rank, int p);
void   vectorToFile(FILE *fout, int s, int n, double* val, int rank, int p);
