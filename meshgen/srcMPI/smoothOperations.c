// ##################################################################
//
// smoothOperations.c
// 
// Associated smoothing operations
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "smoothOperations.h"
#include "globalVariables.h"
#include "meshtype.h"

// ##################################################################
//
// Routines to add the vertices onto the edge of a triangle
//
// posVert = size(3*nVert)
// ##################################################################
void addVerticesOnEdge(const double *posA,      
                       const double *posB,
                             double *posVert,
                       const double *normalA,
                       const double *normalB,
                             double *normalVert,
                       const int     iSurface, 
                       const int     nVert,
                       const double  deltaDist)
{

   // ===============================================================
   // New points on the edge are placed at equal intervals
   // ===============================================================
   int i,j,k;

   k = 0;
   for (i = 0; i < nVert; i++)
   {
   
      posVert[k  ] = posA[0] + (i+1)*deltaDist*(posB[0]-posA[0]);
      posVert[k+1] = posA[1] + (i+1)*deltaDist*(posB[1]-posA[1]);
      posVert[k+2] = posA[2] + (i+1)*deltaDist*(posB[2]-posA[2]);
      
      normalVert[k  ] = normalA[0] + (i+1)*deltaDist*(normalB[0]-normalA[0]);
      normalVert[k+1] = normalA[1] + (i+1)*deltaDist*(normalB[1]-normalA[1]);
      normalVert[k+2] = normalA[2] + (i+1)*deltaDist*(normalB[2]-normalA[2]);
      

      // pass address of the particular vertex only
      if(iSurface == 1)
         moveToBoundary(posA,posB,&posVert[k],&normalVert[k],surfaceType);
      k+=3;

   } // i loop

}

// ################################################################## 
// SmoothGrid
// Inputs
//  - Pointer to the grid
//  - msweep (number of smoothing passes)
// ##################################################################
void smoothGrid(GRID *g, const int msweep)
{
   printf("#meshgen: Lagrangian smoothing ...\n\n");

   int      i,m,j,threeNum;
   int      node1,node2;
   int     *iflag;
   int     *nodeCount;   
   int    **indx;
   double   norm,dummy1[3],dummy2[3]; 
   double  *x1;
   double **weight;
   
   //
   indx       = (int **)   malloc(sizeof(int *)*g->numNodePos);
   iflag      = (int *)    malloc(sizeof(int)*g->numNodePos);
   nodeCount  = (int *)    malloc(sizeof(int)*g->numNodePos);
   x1         = (double *) malloc(sizeof(double)*3*g->numNodePos);

   for (i = 0; i < 3; i++)
      dummy1[i] = dummy2[i] = 0.0;

   // set all iflag and nodeCount for all nodes to 0
   for (i = 0; i < g->numNodePos; i++) iflag[i]=nodeCount[i]=0;

   // additional condition for robin
   if (strcmp(surfaceType,"robin")==0)
   {
      for (i = 0; i < g->numTriNode; i++) iflag[i]=1;
   }

   // loop through all the faces (i.e., edges)
   for (i = 0; i < g->numQuadEdge; i++)
   {
      
      // two node IDs of a given edge
      node1 = g->quadEdge[i][0];
      node2 = g->quadEdge[i][1];
      
      nodeCount[node1]++; 
      nodeCount[node2]++;      

      // Perform swapping operations
      // swap if left cell = -1
      // therefore as a result all -1s are now in [i][3]
      if (g->quadEdge[i][2] == -1 || g->quadEdge[i][2] == -5)
      {
         SWAP(g->quadEdge[i][0],g->quadEdge[i][1]);
         SWAP(g->quadEdge[i][2],g->quadEdge[i][3]);
         SWAP(g->quadEdge[i][4],g->quadEdge[i][5]);
      }
      
      // if the left/right face is -1
      // set iflag to 1
      if (g->quadEdge[i][3] == -1 || g->quadEdge[i][3] == -5)
      {
         iflag[node1] = 1;
         iflag[node2] = 1;
      }



   }
   // indx is a double pointer. So allocate memory for the
   // first level of pointer
   for (i = 0; i < g->numNodePos; i++)
      indx[i]   = (int *)    malloc(sizeof(int)*nodeCount[i]); 
      
   // reset nodeCount array
   for (i = 0; i < g->numNodePos; i++) nodeCount[i]=0;
  
   // The indx pointer points to a pointer of nnodes and each 
   // of pointer of a given node contains the node IDs of the 
   // other end of the edge
   for (i = 0; i < g->numQuadEdge; i++) // loop through all the faces
   {
      node1 = g->quadEdge[i][0];
      node2 = g->quadEdge[i][1];
      
      indx[node1][nodeCount[node1]] = node2;
      nodeCount[node1]++;

      indx[node2][nodeCount[node2]] = node1;
      nodeCount[node2]++;
   }



   // ===============================================================
   // Smoothing procedure
   // Creates a new array x1 (which contain the smoothed
   // values of x and y). The procedure is a simple average
   // where the center node is replaced by the average of 
   // the connecting nodes. If it is a edge node 
   // (i.e., iflag=1), the original positions are retained
   // ===============================================================
   threeNum = 3*g->numNodePos;
   // loop through the number of sweeps
   for(m = 0; m < msweep; m++)
   {          
      // loop through all the nodes    
      for(i = 0; i < g->numNodePos; i++)
      {
         // iflag is 1 only at edges of the boundary
         // iflag is 0 otherwise
         // do not alter original triangular mesh points?!
         if (iflag[i]==0)
         {
            x1[3*i] = x1[3*i+1] = x1[3*i+2] = 0.;
 
            // loop through the number of edges connecting
            // with the current node
            for(j = 0; j < nodeCount[i]; j++)
            {
               x1[3*i]   += g->allNodePos[3*indx[i][j]  ];
               x1[3*i+1] += g->allNodePos[3*indx[i][j]+1];
               x1[3*i+2] += g->allNodePos[3*indx[i][j]+2];
            }
            
            // compute the average x and y position
            x1[3*i  ]/=nodeCount[i];
            x1[3*i+1]/=nodeCount[i];
            x1[3*i+2]/=nodeCount[i];
            
         }
         else
         {
            x1[3*i  ] = g->allNodePos[3*i  ];
            x1[3*i+1] = g->allNodePos[3*i+1];
            x1[3*i+2] = g->allNodePos[3*i+2];
         }
       }
      norm=0.0;

      // Perform smoothing if neccessary
      if (strcmp(surfaceType,"sphere")==0 || 
          strcmp(surfaceType,"robin" )==0)
      {
         for (i = 0; i < g->numNodePos; i++)
            moveToBoundary(dummy1,dummy2,&x1[3*i],
               &g->nodeNormal[3*i],surfaceType);
      }

      // Compute the L2 norm and replace x with x1
      for (i = 0; i < threeNum; i++)
      {

         norm            += (g->allNodePos[i]-x1[i])*(g->allNodePos[i]-x1[i]);
         g->allNodePos[i] = x1[i];
      }
      norm = sqrt(norm)/threeNum;
      //tracef(norm);
   } // m loop (msweep)


   // ===============================================================
   // Free the memory used
   // ===============================================================
  free(nodeCount);
  free(x1);

  for(i = 0; i < g->numNodePos; i++)
    free(indx[i]);

  free(indx);
  free(iflag);
}

// ################################################################## 
// SmoothTriangleGrid
//
// Routines that smoothes only the triangular nodes, midpoints, and
// centroids. The quad subdivision is done after this step.
//
// Inputs
//  - Pointer to the grid
//  - msweep (number of smoothing passes)
// ##################################################################
void smoothTriangleGrid(GRID *g, const int msweep)
{
   printf("#meshgen: Lagrangian smoothing (triangular grid) ...\n\n");

   int      i,m,j,threeNum;
   int      node1,node2;
   int     *iflag;
   int     *nodeCount;   
   int    **indx;
   double   norm,dummy1[3],dummy2[3]; 
   double  *x1;
   double **weight;

   // writemultiarrayINT(g->triQuadEdge,g->numTriQuadEdge,4);
   // traces(stopping);exit(1);
   
   //
   indx       = (int **)   malloc(sizeof(int *)*g->numNodePos);
   iflag      = (int *)    malloc(sizeof(int)*g->numNodePos);
   nodeCount  = (int *)    malloc(sizeof(int)*g->numNodePos);
   x1         = (double *) malloc(sizeof(double)*3*g->numNodePos);

   for (i = 0; i < 3; i++)
      dummy1[i] = dummy2[i] = 0.0;

   // set all iflag and nodeCount for all nodes to 0
   for (i = 0; i < g->numNodePos; i++) iflag[i]=nodeCount[i]=0;

   // additional condition for robin
   if (strcmp(surfaceType,"robin")==0)
   {
      for (i = 0; i < g->numTriNode; i++) iflag[i]=1;
   }
   // loop through all the faces (i.e., edges)
   for (i = 0; i < g->numTriQuadEdge; i++)
   {
      
      // two node IDs of a given edge
      node1 = g->triQuadEdge[i][0];
      node2 = g->triQuadEdge[i][1];

      nodeCount[node1]++; 
      nodeCount[node2]++;      

      // Perform swapping operations
      // swap if left cell = -1
      // therefore as a result all -1s are now in [i][3]
/*
      if (g->quadEdge[i][2] == -1 || g->quadEdge[i][2] == -5)
      {
         SWAP(g->quadEdge[i][0],g->quadEdge[i][1]);
         SWAP(g->quadEdge[i][2],g->quadEdge[i][3]);
         SWAP(g->quadEdge[i][4],g->quadEdge[i][5]);
      }

*/
      // Perform swapping operations
      // swap if left cell = -1
      // therefore as a result all -1s are now in [i][3]
      // ,i.e., to the right side of the edge
      if (g->triQuadEdge[i][2] == -1 || g->triQuadEdge[i][2] == -5)
      {
         SWAP(g->triQuadEdge[i][0],g->triQuadEdge[i][1]);
         SWAP(g->triQuadEdge[i][2],g->triQuadEdge[i][3]);
      }
    
      // if the left/right face is -1
      // set iflag to 1
      if (g->triQuadEdge[i][3] == -1 || g->triQuadEdge[i][3] == -5)
      {
         iflag[node1] = 1;
         iflag[node2] = 1;
      }
   }

   // indx is a double pointer. So allocate memory for the
   // first level of pointer
   for (i = 0; i < g->numNodePos; i++)
      indx[i]   = (int *) malloc(sizeof(int)*nodeCount[i]); 
      
   // reset nodeCount array
   for (i = 0; i < g->numNodePos; i++) nodeCount[i]=0;
  
   // The indx pointer points to a pointer of nnodes and each 
   // of pointer of a given node contains the node IDs of the 
   // other end of the edge
   for (i = 0; i < g->numTriQuadEdge; i++) // loop through all the faces
   {
      node1 = g->triQuadEdge[i][0];
      node2 = g->triQuadEdge[i][1];
      
      indx[node1][nodeCount[node1]] = node2;
      nodeCount[node1]++;

      indx[node2][nodeCount[node2]] = node1;
      nodeCount[node2]++;
   }

   // ===============================================================
   // Smoothing procedure
   // Creates a new array x1 (which contain the smoothed
   // values of x and y). The procedure is a simple average
   // where the center node is replaced by the average of 
   // the connecting nodes. If it is a edge node 
   // (i.e., iflag=1), the original positions are retained
   // ===============================================================
   threeNum = 3*g->numNodePos;
   // loop through the number of sweeps
   for(m = 0; m < msweep; m++)
   {          
      // loop through all the nodes    
      for(i = 0; i < g->numNodePos; i++)
      {
         // iflag is 1 only at edges of the boundary
         // iflag is 0 otherwise
         if (iflag[i]==0 && nodeCount[i]!=0)
         {
            x1[3*i] = x1[3*i+1] = x1[3*i+2] = 0.;
 
            // loop through the number of edges connecting
            // with the current node
            for(j = 0; j < nodeCount[i]; j++)
            {
               x1[3*i  ] += g->allNodePos[3*indx[i][j]  ];
               x1[3*i+1] += g->allNodePos[3*indx[i][j]+1];
               x1[3*i+2] += g->allNodePos[3*indx[i][j]+2];
            }
            
            // compute the average x and y position
            x1[3*i  ]/=nodeCount[i];
            x1[3*i+1]/=nodeCount[i];
            x1[3*i+2]/=nodeCount[i];
         }
         else
         {
            x1[3*i  ] = g->allNodePos[3*i  ];
            x1[3*i+1] = g->allNodePos[3*i+1];
            x1[3*i+2] = g->allNodePos[3*i+2];
         }
      }
      norm=0.0;

      // Perform smoothing if neccessary
      if (strcmp(surfaceType,"sphere")==0 ||
          strcmp(surfaceType,"robin")==0)
      {
         for (i = 0; i < g->numNodePos; i++)
            moveToBoundary(dummy1,dummy2,&x1[3*i],
            &g->nodeNormal[3*i],surfaceType);
      }

      // Compute the L2 norm and replace x with x1
      for (i = 0; i < threeNum; i++)
      {
         norm            += (g->allNodePos[i]-x1[i])*(g->allNodePos[i]-x1[i]);
         g->allNodePos[i] = x1[i];
      }
      norm = sqrt(norm)/threeNum;
      //tracef(norm);
   } // m loop (msweep)


   // ===============================================================
   // Free the memory used
   // ===============================================================
   free(nodeCount);
   free(x1);

   for(i = 0; i < g->numNodePos; i++)
      free(indx[i]);

   free(indx);
   free(iflag);
}

// ##################################################################
//
// smooth the triangle edges
//
// ##################################################################
void smoothTriEdge(GRID *g)
{
   printf("#meshgen: Smoothing triangle edges ...\n");

   int i,i1,i2,j,k1,iSurface;
   double posA[3],posB[3],posA1[3],posA2[3],posA3[3];


   k1 = g->numTriNode;
   // Loop over total number of edges to generate 1/4, 1/2 and 3/4 
   // points along each edge and create basic connectivity info
   for (i = 0; i < g->numTriEdge; i++)
   {
   
      i1 = g->triEdge[i][0]; // one of the end nodes of the dge
      i2 = g->triEdge[i][1]; // other end node of the edge

      //      0.25   0.25    0.25   0.25
      //    |------|------||------|------|
      //    A     A1      A2      A3    B
      for (j = 0; j < 3; j++ )
      {
         posA[j] = g->allNodePos[3*i1 + j];
         posB[j] = g->allNodePos[3*i2 + j];  
      }
   
      // Is this a domain specific condition. Must check.
      if ( (g->triEdge[i][3]==-1) && (abs(posA[1])<0.61) && (abs(posA[0])<1.1) )
         iSurface=1;                  
      else
         iSurface=0;      


      if(iSurface == 0)
      {
         k1 = g->numTriNode + 3*i;

         for (j = 0; j < 3; j++)
         {
            posA1[j] = g->allNodePos[3*k1    +j];
            posA2[j] = g->allNodePos[3*(k1+1)+j];
            posA3[j] = g->allNodePos[3*(k1+2)+j];               
         }
      

         // Averaging the position
         for (j = 0; j < 3; j++)
         {
            posA1[j] = 0.5*(posA[j] + posA2[j]);
            posA3[j] = 0.5*(posB[j] + posA2[j]);
         }

         // save back into allNodePos
         for (j = 0; j < 3; j++)
         {
            g->allNodePos[3*k1    +j] = posA1[j];
            g->allNodePos[3*(k1+2)+j] = posA3[j];               
         }

      }
   } // i loop

}


// ##################################################################
// END OF FILE
// ##################################################################
