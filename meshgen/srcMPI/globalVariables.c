// ##################################################################
//
// globalVariables.c
//
// List of global variables
// ##################################################################
#include "meshtype.h"
#include "globalVariables.h"


int     nGrid;          // total number of grids
int     numSmooth;      // number of smoothing passes on the mesh
int     iStrand;        // are strand grids required
int     pow2;           // 2^(quadlevel)
int     pow4;           // 4^(quadlevel)
char    surfaceType[20];// type of surface to flush the newly formed points to
char    folderName[250];// folder for the subdomain
char    smoothTechnique[20]; // blend or lagrangian
GRID   *g;              // pointer to grid

// for robin fuselage
double *rc; // robin c,d,e,f
double *rd;
double *re;
double *rf;

// MPI related
const int MAX_CONN = 5;
int       totalBoundaryNodeCount;
int       maxNumNodes;
int       maxNodeConnections;
int       totalNumNodePos;
int     **allFwdMap;
int     **allNodeRevMap;
int      *allBoundaryNodeID;
int     **allBoundaryMapID;
int      *allBoundaryNodeCount;
int      *allNumNodePos;
int      *allBoundaryNodeCountCum;
int      *allNumNodePosCum;
int      *subdomainconn;
int      *subdomainOtherID;
double   *allBoundaryNodePos;
// ##################################################################
// enD OF FILE
// ##################################################################
