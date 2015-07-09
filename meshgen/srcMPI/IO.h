// ##################################################################
//
// IO.h
//
// Header file for input output routines
// ##################################################################
#ifndef IO_H
#define IO_H

#include "meshtype.h"

void welcome();

void thanks();

void readInputs(int gridID);

void writeTecplot(int gridID, GRID *g);

void createOutputFolder(int gridID);

void writearrayINT(const int *array, const int n);

void writearrayDUB(const double *array, const int n);

void writemultiarrayINT(const int ** array, const int m, const int n);

void writemultiarrayDUB(const double ** array, const int m, const int n);

#endif
// ##################################################################
// END OF FILE
// ##################################################################
