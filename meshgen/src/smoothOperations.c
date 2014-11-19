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
      
      // pass address of the particular vertex only
      if(iSurface == 1)
         moveToBoundary(posA,posB,&posVert[k],surfaceType);
      k+=3;

   } // i loop

}

// ################################################################## 
// computeNodeWeights
// 
// Computes the weights associated with each node
// ##################################################################
void computeNodeWeights(GRID *g)
{
   printf("#meshgen: Computing node weights ...\n");
   // logic - each node - find triangle loops - 
   //         each triangle - find edges
   //         if node part of an edge - compute associated length
   //         rinse and repeat
   int    i,j,k,i1,i2,nv;
   int    vertID;
   double posA[3],posB[3],sum;

   g->triNodeWeight = (double *) malloc(sizeof(double)*g->numTriNode);
   
   // loop through all nodes
   for (i = 0; i < g->numTriNode; i++)
   {
      posA[0] = g->nodePosTri[3*i  ];
      posA[1] = g->nodePosTri[3*i+1];
      posA[2] = g->nodePosTri[3*i+2];

      sum = 0.0;
      i1 = g->vertLoop->index[i];
      i2 = g->vertLoop->index[i+1];
      nv = i2-i1;
      // loop through the vertices of a particular node
      for (j = i1; j < i2; j++)
      {
         // ID of the vertex
         vertID = g->vertLoop->ID[j];

         // position of vertex
         posB[0] = g->nodePosTri[3*vertID  ];
         posB[1] = g->nodePosTri[3*vertID+1];
         posB[2] = g->nodePosTri[3*vertID+2];

         sum += (posA[0]-posB[0])*(posA[0]-posB[0]) +
                (posA[1]-posB[1])*(posA[1]-posB[1]) +
                (posA[2]-posB[2])*(posA[2]-posB[2]);
      }
      sum = sqrt(sum)/nv;
      g->triNodeWeight[i] = sum;

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
   indx       = (int **)    malloc(sizeof(int *)*g->numNodePos);
   iflag      = (int *)     malloc(sizeof(int)*g->numNodePos);
   nodeCount  = (int *)     malloc(sizeof(int)*g->numNodePos);
   x1         = (double *)  malloc(sizeof(double)*3*g->numNodePos);

   for (i = 0; i < 3; i++)
      dummy1[i] = dummy2[i] = 0.0;

   // set all iflag and nodeCount for all nodes to 0
   for (i = 0; i < g->numNodePos; i++) iflag[i]=nodeCount[i]=0;

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
      if (g->quadEdge[i][2] == -1)
      {
         SWAP(g->quadEdge[i][0],g->quadEdge[i][1]);
         SWAP(g->quadEdge[i][2],g->quadEdge[i][3]);
         SWAP(g->quadEdge[i][4],g->quadEdge[i][5]);
      }
      
      // if the left/right face is -1
      // set iflag to 1
      if (g->quadEdge[i][3] == -1)
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
      //for(i = 400; i < 401; i++)
      {
         // iflag is 1 only at edges of the boundary
         // iflag is 0 otherwise
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
      if (strcmp(surfaceType,"sphere")==0)
      {
         for (i = 0; i < g->numNodePos; i++)
            moveToBoundary(dummy1,dummy2,&x1[3*i],surfaceType);
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
// smooth the loops and not the nodes
//
// ##################################################################
void smoothLoops(GRID *g)
{

   printf("#meshgen: Smoothing along loops ...\n");

   int     i,ii,j,k1,k2,m,kk;
   int     msweep,temp0,temp1,temp2,threeNum;
   int     i1,i2,nel,e1,e2,is,ie;
   int     edge1[2],edge2[2];
   int     *vert4ID,*vert4Index;
   double *x1,sum1,sum2,sum3;
   x1 = (double *) malloc(sizeof(double)*3*g->numNodePos);
   
   threeNum = 3*g->numNodePos;
   for (i = 0; i < threeNum; i++)
      x1[i] = g->allNodePos[i];

   g->vert1Loop = (LOOP *) malloc(sizeof(LOOP));
   g->vert2Loop = (LOOP *) malloc(sizeof(LOOP));
   g->vert3Loop = (LOOP *) malloc(sizeof(LOOP));

   // set the total length of these loops
   g->vert1Loop->totLen = 0.5*(g->quadLoop->totLen);
   g->vert2Loop->totLen = 0.5*(g->quadLoop->totLen);
   g->vert3Loop->totLen = 0.5*(g->quadLoop->totLen);

   g->vert1Loop->ID = (int *) malloc(sizeof(int)*g->vert1Loop->totLen);
   g->vert2Loop->ID = (int *) malloc(sizeof(int)*g->vert2Loop->totLen);
   g->vert3Loop->ID = (int *) malloc(sizeof(int)*g->vert3Loop->totLen);
   vert4ID          = (int *) malloc(sizeof(int)*g->vert3Loop->totLen);

   g->vert1Loop->index = (int *) malloc(sizeof(int)*(1+g->numTriNode));
   g->vert2Loop->index = (int *) malloc(sizeof(int)*(1+g->numTriNode));
   g->vert3Loop->index = (int *) malloc(sizeof(int)*(1+g->numTriNode));
   vert4Index          = (int *) malloc(sizeof(int)*(1+g->numTriNode));

   int iouter  = 1;
   int iinner  = 1;
   int imiddle = 1;
   int iinnest = 1;

   // Initialize
   for (i = 0; i <= g->numTriNode; i++)
   {
      g->vert1Loop->index[i] = g->vert2Loop->index[i] = 0;
      g->vert3Loop->index[i] = vert4Index[i] = 0;
   }

// ==================================================================
// Set up the indices for the different loops
// 4 loops per node
// ==================================================================   
   k1 = 0; k2 = 0; kk = 0;
   // Loop through the number of nodes
   for (i = 0; i < g->numTriNode; i++)
   {
      //ii = g->colourIndex[i];

      // ============================================================
      // INNERMOST LOOP
      // ============================================================
      // indices for start and end of the inner and outer loop
      // for a given triangular node
      i1  = g->iqloop3[i];
      i2  = g->iqloop3[i+1];      
      nel = i2-i1;

      // Loop through the inner
      for (j = i1; j < i2; j++)
      {
         e1 = g->q3loop[j];         

         edge1[0] = g->quadEdge[e1][0];
         edge1[1] = g->quadEdge[e1][1];

         //vert4ID[k2] = edge1[0];
         vert4ID[k2] = edge1[1];
         k2++;
      } // j loop
      vert4Index[kk+1] = vert4Index[kk] + nel;
   
      // ============================================================
      // 3 other loops
      // ============================================================
      // indices for start and end of the inner and outer loop
      // for a given triangular node
      i1  = g->quadLoop->index[2*i  ];
      i2  = g->quadLoop->index[2*i+1];            
      nel = i2-i1;      

      // Loop through the inner
      for (j = i1; j < i2; j++)
      {
         e1 = g->quadLoop->ID[j];
         e2 = g->quadLoop->ID[j+nel];

         edge1[0] = g->quadEdge[e1][0];
         edge1[1] = g->quadEdge[e1][1];

         edge2[0] = g->quadEdge[e2][0];
         edge2[1] = g->quadEdge[e2][1];

         // Note that edge1[1] and edge2[0] are always
         // the same. Therefore the other remaining
         // vertex IDs are those of the two loops
         g->vert1Loop->ID[k1] = edge1[0];
         g->vert2Loop->ID[k1] = edge2[1];
         g->vert3Loop->ID[k1] = edge1[1];
         k1++;
      } // j loop


      g->vert1Loop->index[kk+1] = g->vert1Loop->index[kk] + nel;
      g->vert2Loop->index[kk+1] = g->vert2Loop->index[kk] + nel;
      g->vert3Loop->index[kk+1] = g->vert3Loop->index[kk] + nel;

      
      kk++;


   } // i loop

// ==================================================================
// With the created vert loops - apply laplacian smoothing
// to each of the loops
// Perform the sweeps   
// ==================================================================   
   msweep = 5;
   traces('are you sure that msweep is 5?');
   exit(1);
   for (m = 0; m < msweep; m++)
   {
      for (i = 0; i < g->numTriNode; i++)
      {         
         
         // =========================================================
         // INNERMOST LOOP
         // =========================================================
         i1 = vert4Index[i];
         i2 = vert4Index[i+1];
         nel= i2-i1;
         sum1=sum2=sum3=0.0;nel=0;
         // smooth the positions (same for open or closed loop)
         for (j = i1+1; j < i2-2; j++)
         {  
            temp0 = 3*vert4ID[j];
            temp1 = 3*vert4ID[j-1];
            temp2 = 3*vert4ID[j+1];
            
            if(iinnest==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
               sum1 += g->allNodePos[temp0];
               sum2 += g->allNodePos[temp0+1];
               sum3 += g->allNodePos[temp0+2];
               nel++;
            }

         } // j loop
         
         if( vert4ID[i1] == vert4ID[i2-1])
         {
            temp0 = 3*vert4ID[i1];
            temp1 = 3*vert4ID[i1+1];
            temp2 = 3*vert4ID[i2-1];

            if(iinnest==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
               sum1 += g->allNodePos[temp0];
               sum2 += g->allNodePos[temp0+1];
               sum3 += g->allNodePos[temp0+2];
               nel++;
            }
         
         }
         ii = g->colourIndex[i];
         if(iinnest==1)
         {
            x1[3*ii  ] = sum1/nel;
            x1[3*ii+1] = sum2/nel;
            x1[3*ii+2] = sum3/nel;
         }
         
         // =========================================================
         // OUTER 3 LOOPS
         // =========================================================
         i1  = g->vert1Loop->index[i];
         i2  = g->vert1Loop->index[i+1];
         nel = i2-i1;

         // smooth the positions (same for open or closed loop)
         for (j = i1+1; j < i2-2; j++)
         {  
            
            // OUTER LOOP       
            temp0 = 3*g->vert1Loop->ID[j];
            temp1 = 3*g->vert1Loop->ID[j-1];
            temp2 = 3*g->vert1Loop->ID[j+1];
            
            if(iouter==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }
            // INNER LOOP         
            temp0 = 3*g->vert2Loop->ID[j];
            temp1 = 3*g->vert2Loop->ID[j-1];
            temp2 = 3*g->vert2Loop->ID[j+1];
            
            if(iinner==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }

            // MIDDLE LOOP
            temp0 = 3*g->vert3Loop->ID[j];
            temp1 = 3*g->vert3Loop->ID[j-1];
            temp2 = 3*g->vert3Loop->ID[j+1];
            
            if(imiddle==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }
         } // j loop
         
         // check if closed loop or not. Done by checking the
         // vert ID of the beginning and end of the loop. Can use
         // wither vert1Loop or vert2Loop      
         
         if( g->vert1Loop->ID[i1] == g->vert1Loop->ID[i2-1])
         {
            temp0 = 3*g->vert1Loop->ID[i1];
            temp1 = 3*g->vert1Loop->ID[i1+1];
            temp2 = 3*g->vert1Loop->ID[i2-1];

            if(iouter==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }

            // INNER LOOP
            temp0 = 3*g->vert2Loop->ID[i1];
            temp1 = 3*g->vert2Loop->ID[i1+1];
            temp2 = 3*g->vert2Loop->ID[i2-1];
            
            if(iinner==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }

            // MIDDLE LOOP
            temp0 = 3*g->vert3Loop->ID[i1];
            temp1 = 3*g->vert3Loop->ID[i1+1];
            temp2 = 3*g->vert3Loop->ID[i2-1];
            
            if(imiddle==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }
         }

      } // i loop

// ==================================================================         
// POSITION UPDATE
// ==================================================================            
      for (i = 0; i < g->numNodePos; i++)
      {         
         //i = g->colourIndex[ii];
         g->allNodePos[3*i  ] = g->iBoundary[i]*g->allNodePos[3*i  ]+x1[3*i  ]*(1-g->iBoundary[i]);
         g->allNodePos[3*i+1] = g->iBoundary[i]*g->allNodePos[3*i+1]+x1[3*i+1]*(1-g->iBoundary[i]);
         g->allNodePos[3*i+2] = g->iBoundary[i]*g->allNodePos[3*i+2]+x1[3*i+2]*(1-g->iBoundary[i]);
      }
      
   } // m loop

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
