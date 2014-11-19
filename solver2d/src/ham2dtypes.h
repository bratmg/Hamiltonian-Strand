// ##################################################################
//
// ham2dtypes.h
//
// Code 
// ##################################################################
typedef struct faceMat
{
  double lmat[4][4];
  double rmat[4][4];
} faceMat;

// ==================================================================
typedef struct GRID
{
  int    ncells;  // total number of quadrilatera
  int    nnodes;  // total number of nodes
  int    nfaces;  // total number of faces
  int    nbfaces; // total number of boundary faces
  
  int    ncolors; // total number of colors for chains
  int    nchains; // total number of chains
  int    nchainFaces; // total number of faces in the chainConn list
                   // this is different from nfaces because
                   // some faces are repeated (as they form loops)
                   // in chainConn

  double *x;                 // coordinates of nodes
  int    *conn;              // connectivity of nodes
  int    *chainsPerColor;    // number of chains per color
  int    *faceStartPerChain; // number of faces per chain
  int    *chainConn;         // face connectivity for each chain
                             // stacked together based on numbers in facesPerChain[]

  int    *faces;  // list of faces
  int    *neig;   // list of neighbors
  int    *bfaces; // list of boundary faces
  double *vol; // cell volume
  int    order;   // solution order on this grid
  int    timeacc; // time accurate or not
  int    timeInteg; // timeIntegration scheme
  double CFL;    // nominal CFL number
  int    visc;   // viscous 
  int    test;
  int    msweep;
  // 
  // work arrays for processing fluxes
  // 
  int    nmaxchain;
  double **f;
  double **fv;
  double **ql;
  double **qr;
  double **dqr;
  double **dql;
  double **df;
  double **f2;
  double ***A;
  double ***B;
  double ***C;
  double **F;
  double **Q;
  int    *cindx;
  int    *ctype;
  faceMat *ff;
  int    *c2f;
  int    *c2chain;

  int bf1[5000];
  int bf2[5000];
  int bf3[5000];
  int bf4[5000];
  int m1,m2,m3,m4;

  int bf1_neigf[5000][3];
  int bf2_neigf[5000][3];
  int bf3_neigf[5000][3];
  int bf4_neigf[5000][3];


} GRID;

// ==================================================================
typedef struct SOLN
{
  double *q;     // q -variables [rho rho*u rho*v e]
  double *dq;    // \delta q -variables [rho rho*u rho*v e]
  double *ddq;   // \delta\delta q - variables [rho rho*u rho*v e]
  double *ddqb;  // \delta\delta q - variables [rho rho*u rho*v e]
  double *ddqf;  // \delta\delta q - variables [rho rho*u rho*v e]
  double *r;     // solution residual at cell centroids
  double *r0;    // solution residual at cell centroids
  double *sigma; // line integral of spectral radius per cell
  double ***D;   // diagonal matrix
  double *dtac;  // time scaling - nominally equaly to dt/dv but can be 
                 // combination of pseudo time step as well
  int *itag;
  double l2norm;
  double uinf,vinf,einf;
  double mach,alpha,rey;
  double gm1,c2b,rgas,pr,prtr;
  double cl,cd,cy,cx;
  double *qt;
  double *qtt;
  double cflnum;
  double res0;

} SOLN;
  
// ==================================================================
# define tracef(x) printf("#ham2d:\t"#x" =%.16e\n",x);
# define trace(x)  printf("#ham2d:\t"#x" =%d\n",x);
# define max(x,y)  (x) > (y) ? (x) : (y)
# define min(x,y)  (x) < (y) ? (x) : (y)
# define swap(x,y)  {(x)=(x)+(y);(y)=(x)-(y);(x)=(x)-(y);}
// ##################################################################
// END OF FILE
// ##################################################################
