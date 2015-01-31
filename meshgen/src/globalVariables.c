// ##################################################################
//
// globalVariables.c
//
// List of global variables
// ##################################################################
#include "meshtype.h"

int     nGrid;          // total number of grids
int     numSmooth;      // number of smoothing passes on the mesh
int     iStrand;        // are strand grids required
int     pow2;           // 2^(quadlevel)
int     pow4;           // 4^(quadlevel)
char    surfaceType[20];// type of surface to flush the newly formed points to
char    folderName[250];// folder for the subdomain
GRID   *g;              // pointer to grid

// for robin fuselage
double *rc; // robin c,d,e,f
double *rd;
double *re;
double *rf;

// ##################################################################
// enD OF FILE
// ##################################################################
