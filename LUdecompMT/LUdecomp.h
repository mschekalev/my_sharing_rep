struct Params {
 	int   size;
	int   scale;
    int   threads;
	char  verbose;
	char* fin_name;
	char* fout_name;
	char* formula;
};

typedef struct {
	int     n;
	double* mat;
	double* val;
	double* sol;
    double  time;
    int     fl;
	int     rank;
	int     threads;
} Args;


double currentTimeSec(void);
double currentTimeSecThread(void);
double normFormula(int n, double* a);
double normFile(int n, double* a, double* b);
void*  Solution(void *pa);
void   multiplicationMatVec(int n, double* mat, double* v, double* res);
void   matrixToFile(FILE* fout, int s, int n, double* mat);
void   synchronize(int total_threads);
int    systemSolution(int n, double* mat, double* val, double* sol, int rank, int threads, int *fl);
int    matrixFromFile(FILE* fin, int n, double* mat, double* mat1, double* val);
int    matrixFromFormula(int n, char* formula, double* mat, double* val);
