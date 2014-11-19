// ##################################################################
//
// smoothOperations.h
// 
// Header files for all smoothing operations
// ##################################################################
#ifndef SMOOTH_OPERATIONS_H
#define SMOOTH_OPERATIONS_H

#include "meshtype.h"

void addVerticesOnEdge(const double *posA, const double *posB,
   double *posVert, const int iSurface,
   const int nVert, const double deltaDist);

void computeNodeWeights(GRID *g);

void smoothGrid(GRID *g, const int msweep);

void smoothLoops(GRID *g);

void smoothTriEdge(GRID *g);

#endif
// ##################################################################
// END OF FILE
// ##################################################################
