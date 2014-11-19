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
GRID   *g;              // pointer to grid

// ##################################################################
// enD OF FILE
// ##################################################################
