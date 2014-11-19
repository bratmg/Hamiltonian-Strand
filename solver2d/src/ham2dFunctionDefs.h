//
// definitions for now
//
#define rinf 1.0
#define gamm 1.4
#define pinf 1./gamm
#define deg2rad M_PI/180
#define TWOTHIRD 0.66666666666666667
#define BDF1 1
#define BDF2 2
#define NVAR 4
#define TOL 1e-11
//#define CFL 20.0
#define CFLMAX 10.00
#define RAMPSTART 300000
#define RAMPEND  400000
#define CHK 16350
//
// function prototypes
//
void readGrid(GRID *g);

void preprocess(GRID *g);

void periodic_bc(GRID *g);

void apply_periodic(GRID *g, int f1, int f2, int m);

void apply_periodic_LHS(GRID *g, int f1, int f2, int m);

void initflow(GRID *g,SOLN *s);

void outputSolution(GRID *g,SOLN *s,int nn);

void computeForce(GRID *g,SOLN *s);

void updateSoln(double *qsrc,double *qdest, double *res,double *sigma,
		double dt,double CFL,double coef,int ncells);

void stepSolution(char *stepType,GRID *g,SOLN *s,double dt,double *l2norm);

void muscld(double **f,
	    double **ql,
	    double **qr,
	    double **f2,
	    int is,
	    int ie,
	    double th,
	    double qt,
	    double eps,
	    int imax,
	    int nq);

void muscld_deriv(double **f,
	    double **ql,
	    double **qr,
	    double **f2,
	    double **dq,
	    int is,
	    int ie,
	    double th,
	    double qt,
	    double eps,
	    int imax,
	    int nq);

void roeflx(double *specRadius,double flux[NVAR],
	    double leftState[NVAR],double rightState[NVAR],
	    double faceVel,double ds[2],double gm1);

void smoothGrid(GRID *g, int msweep);

extern void jac_roe_(double *,double *, double *,double **,double **,double *,
	             int *);

extern void jac_roe2d_(double *,double *, double *,double **,double **,double *,
	             int *);

extern void flux_roe_(double *,double *,double *,double *,double *,double *);

extern void flux_roe2d_(double *,double *,double *,double *,double *,double *);

extern void wallflux_(double *,double *,double *,double *,double *);

extern void wallfluxjacobian_(double *,double *,double **,double **,double *);

extern void flux_visc_2d_(double *,double *,double *,double *,double *,
			  double *, double *, double *, double *, double *, double *);

void blockTridag(double ***a,double ***b, double ***c,double **f, int N,int nq) ;

void blockTridagPeriodic(double ***a,double ***b, double ***c,double **f, int N,int nq);

void blockTridag4(double ***a,double ***b, double ***c,double **f, int N,int nq) ;

void blockTridagPeriodic4(double ***a,double ***b, double ***c,double **f, int N,int nq);

void computeRHS(GRID *g,SOLN *s,double *l2rho);

void computeRHSk(GRID *g,SOLN *s,double *l2rho);

void computeRHSkv(GRID *g,SOLN *s,double *l2rho);

void computeRHSkcheck(GRID *g,SOLN *s,double *l2rho);

void computeRHS1(GRID *g,SOLN *s,double *l2rho);

void computeRHS2(GRID *g,SOLN *s,double *l2rho);

void computeRHS3(GRID *g,SOLN *s,double *l2rho,double eps2);

void ADI(GRID *g,SOLN *s,double cflnum,double dt);

void DADI(GRID *g,SOLN *s,double cflnum,double dt);

void DADIk(GRID *g,SOLN *s,double cflnum,double dt);

void findDiagonals(GRID *g,SOLN *s,double cflnum,double dt);

void DADIk1(GRID *g,SOLN *s,double cflnum,double dt);

void computeLinearRHS(GRID *g,SOLN *s,double cflnum,double *l2rho);

void computeLinearRHScheck(GRID *g,SOLN *s,double cflnum,double *l2rho);

void computeLHS1(GRID *g,SOLN *s,double dt);

void computeLHS2(GRID *g,SOLN *s,double dt);

void computeTimeScaling(GRID *g, SOLN *s,double cflnum,double dt,int istep);

void addTemporalSource(GRID *g, SOLN *s,double cflnum, double dt, int istep);

void updateTime(GRID *g, SOLN *s);

void invertMat5(double A[5][5],
		double f[5],
		double x[5]);

void invertMat4(double A[4][4],
		double f[4],
		double x[4]);

void axb(double A[4][4],double *x,double *x0,double *b,double fac,int N);

void axb1(double A[4][4],double *x,double *b,double fac,int N);

void gaussSeidel(GRID *g,SOLN *s,double cflnum,double dt);

void lineGaussSeidel(GRID *g,SOLN *s,double cflnum,double dt);

void lineGaussSeidel1(GRID *g,SOLN *s,double cflnum,double dt);

void outputdq(GRID *g,SOLN *s);
