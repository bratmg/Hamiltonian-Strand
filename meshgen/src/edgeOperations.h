// ##################################################################
//
// edgeOperations.h
// 
// Header files for all edge operations
// ##################################################################
#ifndef EDGE_OPERATIONS_H
#define EDGE_OPERATIONS_H

#include "meshtype.h"

void findTriangleEdges(GRID *g);

void createVerticesOnEdge(GRID *g);

void createInteriorVertices(GRID *g);

void subdomainEdgeBlanking(GRID *g);

#endif
// ##################################################################
// END OF FILE
// ##################################################################
