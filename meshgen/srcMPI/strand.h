// ##################################################################
//
// strand.h
// 
// Header files that contain the routines for strand grids
// ##################################################################
#ifndef STRAND_H
#define STRAND_H

#include "meshtype.h"

void createStrands(int gridID, GRID *g);

void createStrandTemplate(int gridID, GRID *g);

void computeNormals(int gridID, GRID *g);

void computeNormals_MPI(int gridID);

void averageNormals(int gridID, GRID *g);

void averageNormals_MPI(int gridID);

void updateVariables(int gridID, GRID *g);

#endif
// ##################################################################
// END OF FILE
// ##################################################################
