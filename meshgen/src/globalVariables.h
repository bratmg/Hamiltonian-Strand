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
extern int     iHybrid;
extern int     iHybPeriodic; 
extern int     nHGrid; 
extern int     numSmooth;
extern int     iStrand;
extern int     qLevel;
extern int     pow2;
extern int     pow4;
extern char    surfaceType[20];
extern char    folderName[250];
extern char    smoothTechnique[20];
extern GRID   *g;
extern HGRID  *h;


// for robin fuselage
extern double *rc; // robin c,d,e,f
extern double *rd;
extern double *re;
extern double *rf;

#endif
// ##################################################################
// END OF FILE
// ##################################################################
