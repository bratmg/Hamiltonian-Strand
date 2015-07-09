// ##################################################################
//
// hybridOperations.h
// 
// Routines for patching the body-fitted structured and outer
// unstructured routines
// ##################################################################
#ifndef HYBRIDOPERATIONS_H
#define HYBRIDOPERATIONS_H

#include "meshtype.h"
#include <stdlib.h>
#include <stdio.h>

void meshPatching(HGRID *h);

void findValue(const int *array, const int n, const int val, int *nodeID);

#endif
// ##################################################################
// END OF FILE
// ##################################################################