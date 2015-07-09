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
extern char    smoothTechnique[20]; // blend or lagrangian
extern GRID   *g; // pointer to grid


// for robin fuselage
extern double *rc; // robin c,d,e,f
extern double *rd;
extern double *re;
extern double *rf;

// MPI related
extern const int MAX_CONN;
extern int      totalBoundaryNodeCount;
extern int      maxNumNodes;
extern int      maxNodeConnections;
extern int      totalNumNodePos;
extern int    **allFwdMap;
extern int    **allNodeRevMap;
extern int     *allBoundaryNodeID;
extern int    **allBoundaryMapID;
extern int     *allBoundaryNodeCount;
extern int     *allNumNodePos;
extern int     *allBoundaryNodeCountCum;
extern int     *allNumNodePosCum;
extern int     *subdomainconn;
extern int     *subdomainOtherID;
extern double  *allBoundaryNodePos;

#endif
// ##################################################################
// END OF FILE
// ##################################################################
