// ##################################################################
//
// smoothGrid.c
//
// Code to smoothen the grid obtained using Hamilton loop generation
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

// ==================================================================  
// SmoothGrid
// Inputs
//  - Pointer to the grid
//  - msweep (number of smoothing passes)
// ==================================================================
void smoothGrid(GRID *g, int msweep)
{
  int    **indx;
  int    *nodeCount;
  int    *iflag;
  int    i,m,j;
  int    node1,node2;
  double *x1;
  double norm; 
  //
  indx      = (int **)   malloc(sizeof(int *)*g->nnodes);
  iflag     = (int *)    malloc(sizeof(int)*g->nnodes);
  nodeCount = (int *)    malloc(sizeof(int)*g->nnodes);
  x1        = (double *) malloc(sizeof(double)*2*g->nnodes);

  // set all iflag and nodeCount for all nodes to 0
  for (i=0;i<g->nnodes;i++) iflag[i]=nodeCount[i]=0;

  // loop through all the faces (i.e., edges)
  for (i=0;i<g->nfaces;i++)
  {
      // two node IDs of a given edge
      node1 = g->faces[6*i];
      node2 = g->faces[6*i+1];

      nodeCount[node1]++; 
      nodeCount[node2]++;      
      
      // if the left/right face is -1
      // set iflag to 1
      if (g->faces[6*i+4] == -1) 
      {
         iflag[node1] = 1;
         iflag[node2] = 1;
      }
  }
  
  // indx is a double pointer. So allocate memory for the
  // first level of pointer
  for (i=0;i<g->nnodes;i++)
    indx[i] = (int *) malloc(sizeof(int)*nodeCount[i]);
  
  // reset nodeCount array
  for (i=0;i<g->nnodes;i++) nodeCount[i]=0;
  
  // The indx pointer points to a pointer of nnodes and each 
  // of pointer of a given node contains the node IDs of the 
  // other end of the edge
  for (i=0;i<g->nfaces;i++) // loop through all the faces
  {
      node1 = g->faces[6*i];
      node2 = g->faces[6*i+1];
      
      indx[node1][nodeCount[node1]] = node2;
      nodeCount[node1]++;

      indx[node2][nodeCount[node2]] = node1;
      nodeCount[node2]++;
  }

// ==================================================================
// Smoothing procedure
// Creates a new array x1 (which contain the smoothed
// values of x and y). The procedure is a simple average
// where the center node is replaced by the average of 
// the connecting nodes. If it is a edge node 
// (i.e., iflag=1), the original positions are retained
// ==================================================================
  // loop through the number of sweeps
  for(m=0;m<msweep;m++)
  {      
    
    // loop through all the nodes    
    for(i=0;i<g->nnodes;i++)
    {
      // iflag is 1 only at edges of the boundary
      // iflag is 0 otherwise
      if (iflag[i]==0)
      {
        x1[2*i] = x1[2*i+1] = 0;

        // loop through the number of edges connecting
        // with the current node
        for(j=0;j<nodeCount[i];j++)
        {
           x1[2*i]   += g->x[2*indx[i][j]];
           x1[2*i+1] += g->x[2*indx[i][j]+1];
        }

        // compute the average x and y position
        x1[2*i]/=nodeCount[i];
        x1[2*i+1]/=nodeCount[i];
      }
      else
      {
        x1[2*i]   = g->x[2*i];
        x1[2*i+1] = g->x[2*i+1];
      }
    }
    norm=0.0;

    // Compute the L2 norm and replace x with x1
    for(i=0;i<2*g->nnodes;i++)
    {
      norm   += (g->x[i]-x1[i])*(g->x[i]-x1[i]);
      g->x[i] = x1[i];
    }
    norm = sqrt(norm)/2/g->nnodes;
    //tracef(norm);
  }

// ==================================================================
// Free the memory used
// ==================================================================
  free(nodeCount);
  free(x1);

  for(i=0;i<g->nnodes;i++)
    free(indx[i]);

  free(indx);
  free(iflag);
}
// ##################################################################
// END OF FILE
// ##################################################################