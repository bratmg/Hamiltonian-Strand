// ##################################################################
//
// globalVariables.c
//
// List of global variables
// ##################################################################
#include "meshtype.h"

int     nGrid;          // total number of grids
int     iHybrid;        // if hybrid grid exists
int     iHybPeriodic;   // yes/no periodic BC for structured grid
int     nHGrid;         // number of hybrid grids
int     numSmooth;      // number of smoothing passes on the mesh
int     iStrand;        // are strand grids required
int     qLevel;         // quad level of subdivision
int     pow2;           // 2^(quadlevel)
int     pow4;           // 4^(quadlevel)
char    surfaceType[20];// type of surface to flush the newly formed points to
char    folderName[250];// folder for the subdomain
char    smoothTechnique[20]; // blend or lagrangian
GRID   *g;              // pointer to grid
HGRID  *h;              // pointer to hybrid grid

// for robin fuselage
double *rc; // robin c,d,e,f
double *rd;
double *re;
double *rf;

// ##################################################################
// enD OF FILE
// ##################################################################
