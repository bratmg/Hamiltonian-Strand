// ##################################################################
//
// meshQuality.h
//
// Header file for routines related to checking mesh quality
// ##################################################################
#ifndef MESH_QUALITY_H
#define MESH_QUALITY_H

#include "meshtype.h"

void meshQuality(int gridID, GRID *g);

void checkQuality(int gridID, GRID *g);

void outputMeshQuality(GRID *g);

#endif
// ##################################################################
// END OF FILE
// ##################################################################
