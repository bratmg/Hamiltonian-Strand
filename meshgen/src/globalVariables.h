// ##################################################################
//
// globalVariables.h
//
// List of global variables
// ##################################################################
#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H

#include "meshtype.h"

extern int     nGrid;
extern int     numSmooth;
extern int     iStrand;
extern int     pow2;
extern int     pow4;
extern char    surfaceType[20];
extern char    folderName[250];
extern GRID   *g; // pointer to grid

// for robin fuselage
extern double *rc; // robin c,d,e,f
extern double *rd;
extern double *re;
extern double *rf;

#endif
// ##################################################################
// END OF FILE
// ##################################################################
