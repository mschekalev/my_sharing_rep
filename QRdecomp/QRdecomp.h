struct Params {
 	int    size;
	int    scale;
	char   verbose;
	char*  fin_name;
	char*  fout_name;
	char*  formula;
    double eps;
};


//double currentTimeSec(void);
double normTrace(int n, double* mat, double* sol);
double normMatrix(int n, double* mat, double* sol);
double qrIteration(int n, int acs, double* mat, double* vcos, double* vsin, double eps);
void   matrixToFile(FILE* fout, int s, int n, double* mat);
void   toAlmostTriangular(int n, double* mat, double eps);
void   qrForAlmostTriangular(int n, int acs, double* mat, double* vcos, double* vsin, double eps);
void   rqMultiplication(int n, int acs, double* mat, double* vcos, double* vsin);
void   qrAlgorithm(int n, double* mat, double* vcos, double* vsin, double* sol, double eps);
void   vectorToFile(FILE* fout, int s, int n, double* vec);
void   matrixShift(int n, int acs, double* mat, double sh);
int    matrixFromFile(FILE* fin, int n, double* mat, double* mat1);
int    matrixFromFormula(int n, char* formula, double* mat, double* mat1);
