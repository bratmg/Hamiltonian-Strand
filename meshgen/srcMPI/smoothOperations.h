// ##################################################################
//
// smoothOperations.h
// 
// Header files for all smoothing operations
// ##################################################################
#ifndef SMOOTH_OPERATIONS_H
#define SMOOTH_OPERATIONS_H

#include "meshtype.h"

void addVerticesOnEdge(const double *posA, 
                       const double *posB,
                             double *posVert, 
                       const double *normalA,
                       const double *normalB,
                             double *normalVert,
                       const    int  iSurface,
                       const    int  nVert, 
                       const double  deltaDist);

void smoothGrid(GRID *g, const int msweep);

void smoothTriangleGrid(GRID *g, const int msweep);

void smoothTriEdge(GRID *g);

#endif
// ##################################################################
// END OF FILE
// ##################################################################
