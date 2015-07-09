// #################################################################
//
// communication.h
//
// Routines that contain MPI based communication routines
// #################################################################
#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <mpi.h> 
#include "meshtype.h"
#include "globalVariables.h"

void boundaryNodeConnection_MPI(int gridID);

#endif
// #################################################################
// END OF FILE
// #################################################################
