// ##################################################################
//
// flushToSurface.h
// 
// Header files that contain the routines to ensure that the
// newly created points are along the surface of the boundary
// ##################################################################
#ifndef FLUSH_TO_SURFACE_H
#define FLUSH_TO_SURFACE_H

void moveToBoundary(double *ptLeft, double *ptRight, 
                    double * pt, char * surfaceID);


#endif
// ##################################################################
// END OF FILE
// ##################################################################
