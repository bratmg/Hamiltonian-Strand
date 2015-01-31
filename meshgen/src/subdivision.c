// ##################################################################
//
// subdivision.c
// 
// Associated edge operations
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "subdivision.h"
#include "globalVariables.h"
#include "meshtype.h"
#include "flushToSurface.h"
#include "smoothOperations.h"

#define N_EDGE 4
#define N_QEDGE 6
#define ONE_THIRD 0.3333333333333333

// ##################################################################
//
// createVerticesOnEdge
//
// Function to create vertices on the triangle edges to form the 
// basis for the different quad cells
//
// ##################################################################
void recreateVerticesOnEdge(GRID *g)
{

   printf("#meshgen: Recreating vertices on triangular edges ...\n");

   int     i,i1,i2,j,k1,k2,ktemp;
   int     iSurface,nVert,nEdge;
   int    *vertID;
   double  posA[3],posB[3],posMid[3],deltaDist;
   double  normalA[3],normalB[3],normalMid[3];
   double *posVert, *normalVert;
   int     halfEdgePerSide = 0.5*g->nEdgePerSide;
   int     halfVertm1      = 0.5*(g->nVertPerSide-1);

   k1 = g->numTriNode;

   // allocate posVert   
   deltaDist  = (double) 2./g->nEdgePerSide;
   posVert    = (double *) malloc(sizeof(double)*3*g->nVertPerSide);
   normalVert = (double *) malloc(sizeof(double)*3*g->nVertPerSide);

   // Loop over total number of edges to generate 1/4, 1/2 and 3/4 
   // points along each edge and create basic connectivity info
   for (i = 0; i < g->numTriEdge; i++)
   {
      i1 = g->triEdge[i][0]; // one of the end nodes of the dge
      i2 = g->triEdge[i][1]; // other end node of the edge

      // index Id for points on edges
      ktemp = k1;

      //      0.25   0.25    0.25   0.25
      //    |------|------||------|------|
      //    A     A1      A2      A3    B
      for (j = 0; j < 3; j++ )
      {
         posA[j]   = g->allNodePos[3*i1 + j];
         posB[j]   = g->allNodePos[3*i2 + j];
         posMid[j] = g->allNodePos[3*(ktemp+halfVertm1) + j];  

         normalA[j]   = g->nodeNormal[3*i1 + j];
         normalB[j]   = g->nodeNormal[3*i2 + j];
         normalMid[j] = g->nodeNormal[3*(ktemp+halfVertm1) + j];  

         if(posMid[j] != posMid[j]) {traces(nan);exit(1);}
      }

      // Is this a domain specific condition. Must check.
      if ( strcmp(surfaceType,"naca")==0
            && (g->triEdge[i][3]==-1) 
            && (abs(posA[1])<0.61) && (abs(posA[0])<1.1) )
         iSurface = 1;                  
      else if (strcmp(surfaceType,"sphere")==0 ||
               strcmp(surfaceType,"robin") ==0 )
         iSurface = 1;
      else
         iSurface = 0;      

      // need to split this into A-D and then D-B
      // as D is altered based on the smoothing routines
      // add new vertices to the edge of a triangle      
      addVerticesOnEdge(posA,posMid,posVert,
                        normalA,normalMid,normalVert,
                        iSurface,halfVertm1,deltaDist);

      // append the position of the newly formed vertices onto
      // the edge array
      for (j = 0; j < halfVertm1; j++)
      {
         g->allNodePos[3*ktemp    ] = posVert[3*j    ];
         g->allNodePos[3*ktemp + 1] = posVert[3*j + 1];
         g->allNodePos[3*ktemp + 2] = posVert[3*j + 2];

         g->nodeNormal[3*ktemp    ] = normalVert[3*j    ];
         g->nodeNormal[3*ktemp + 1] = normalVert[3*j + 1];
         g->nodeNormal[3*ktemp + 2] = normalVert[3*j + 2];

         ktemp++;
      } // j loop
  
      ktemp++; // to skip mid point

      // add new vertices to the edge of a triangle
      addVerticesOnEdge(posMid,posB,posVert,
                        normalMid,normalB,normalVert,
                        iSurface,halfVertm1,deltaDist);

      for (j = 0; j < halfVertm1; j++)
      {
         g->allNodePos[3*ktemp    ] = posVert[3*j    ];
         g->allNodePos[3*ktemp + 1] = posVert[3*j + 1];
         g->allNodePos[3*ktemp + 2] = posVert[3*j + 2];

         g->nodeNormal[3*ktemp    ] = normalVert[3*j    ];
         g->nodeNormal[3*ktemp + 1] = normalVert[3*j + 1];
         g->nodeNormal[3*ktemp + 2] = normalVert[3*j + 2];

         ktemp++;
      } // j loop

      k1 += g->nVertPerSide; 

   } // i loop

   free(posVert);

}

// ##################################################################
//
// createInteriorVertices
//
// Go through the triangles and create the quad cells
//
// ##################################################################
void recreateInteriorVertices(GRID * g)
{

   printf("#meshgen: Recreating interior vertices ...\n");

   int      i,j,jj,n,nn,is1,is2,rowID,colID,counter;
   int      index,halfEdgePerSide,halfVertm1,itemp;
   int      k,kk,k1,k2,k3,k4,k5,k6,ktemp,ktemp3,ktemp4,ktemp5,ktemp6,koffset;
   int      iA,iB,iC,iD,iE,iF,iO;
   int      index1,index2,index3;
   int      edgeID[3],edge[4],dirAB,dirBC,dirCA;
   int     *vertID, iSurface;
   int      nodeMidDO,nodeMidEO,nodeMidFO;
   double   aa,bb,cc;
   double   A[3],B[3],C[3],O[3];
   double   deltaDist,deltaDistTemp;
   double  *boundaryPts, *intPts, *posVert, *normalVert;
   double   dummy1[3],dummy2[3];

   for (i = 0; i < 3; i++)
      dummy1[i] = dummy2[i] = 0.0;

   // ===============================================================
   // Initialization
   // ===============================================================
   halfEdgePerSide = 0.5*g->nEdgePerSide;
   halfVertm1      = 0.5*(g->nVertPerSide-1);
   int    edgePerQuad = 2*halfEdgePerSide*halfVertm1;
   int    halfhalfVertm1  = 0.5*(halfVertm1-1);
   int    halfVertp2 = halfVertm1 + 2;
   int    quadPtID[halfVertp2][halfVertp2];
   double quadPtPos[3*halfVertp2][3*halfVertp2];

   k1 = g->numTriNode + g->nVertPerSide*g->numTriEdge; // for center
   k2 = 0; // for edges
   k3 = g->nEdgePerSide*g->numTriEdge;
   k4 = k1 + g->numTriangle; // int edge points
   k5 = k4 + 3*halfVertm1*g->numTriangle; // quad int points
   k6 = k5; // for additional subdivision at triangular level

   aa = 0.5;
   bb = 0.5;
   cc = 0.5;

   // allocate posVert   
   deltaDist   = 1./halfEdgePerSide;
   posVert     = (double *) malloc(sizeof(double)*3*halfVertm1);
   normalVert  = (double *) malloc(sizeof(double)*3*halfVertm1);

   int     threenum = 3*g->numNodePos;
   double  weight;
   double *tempPos;

   tempPos = (double *) malloc(sizeof(double)*threenum);

   // initialization
   for (i = 0; i < threenum; i++)
      tempPos[i] = g->allNodePos[i];

   // ===============================================================
   //
   // Loop over all the triangles and build the cell connectivity
   // and edge information
   //
   // ===============================================================
   for (i = 0; i < g->numTriangle; i++)
   {
      iA = g->triConn[3*i  ];
      iB = g->triConn[3*i+1];
      iC = g->triConn[3*i+2];

      for (j = 0; j < 3; j++)
      {
         A[j] = g->allNodePos[3*iA + j]; // [x,y,z] coordinate of point A
         B[j] = g->allNodePos[3*iB + j]; // [x,y,z] coordinate of point B
         C[j] = g->allNodePos[3*iC + j]; // [x,y,z] coordinate of point C
      }
      
		// list of the edge ID for a given triangle
		edgeID[0] = g->edge2triList[3*i  ];
		edgeID[1] = g->edge2triList[3*i+1];
		edgeID[2] = g->edge2triList[3*i+2];    

		// Find the coordinate and indices of points on triangle edges
      // Loop over each of the three edges of a triangle
		for (j = 0; j < 3; j++)
		{
         // list columns 1--4 of a particular edge
         edge[0] = g->triEdge[edgeID[j]][0];
         edge[1] = g->triEdge[edgeID[j]][1];         

         // edge index in the allNodePos array
         is1     = g->numTriNode + g->nVertPerSide*(edgeID[j]);
         is2     = g->nEdgePerSide*(edgeID[j]); 

         // Obtain the node ID's of mid points and the direction of 
         // the edge (based on conventions previously defined)
         // =========================================================
         //  If the edge has the 1st node as node A
         // =========================================================
         if (edge[0] == iA)
         {
            // if the edge has the 2nd node as node (i.e, A-B)
            if (edge[1] == iB)
            {
               dirAB = 1;
               iD    = is1 + halfVertm1;
            }
            else if (edge[1] == iC)
            {
               dirCA = -1;
               iF    = is1 + halfVertm1;
            }
         } 

         // =========================================================
         //  If the edge has the 1st node as node B
         // =========================================================
         if (edge[0] == iB)
         {
            // if the edge has the 2nd node as node (i.e, A-B)
            if (edge[1] == iA)
            {
               dirAB = -1;
               iD    = is1 + halfVertm1;
            }
            else if (edge[1] == iC)
            {
               dirBC = 1;
               iE    = is1 + halfVertm1;
            }
         } 

         // =========================================================
         //  If the edge has the 1st node as node C
         // =========================================================
         if (edge[0] == iC)
         {
            // if the edge has the 2nd node as node (i.e, C-A)
            if (edge[1] == iA)
            {
               dirCA = 1;
               iF    = is1 + halfVertm1;  
            }
            else if (edge[1] == iB) // edge is C--B
            {
               dirBC = -1;
               iE    = is1 + halfVertm1;
            }
         }

      } // j loop (edgeID)

      // ============================================================
      // Add the center point
      // ============================================================
      ktemp = k1;
      iO    = ktemp;
      for (k = 0; k < 3; k++)
      {
          O[k] = ONE_THIRD*(A[k] + B[k] + C[k]);
          //g->allNodePos[3*iO + k] = O[k];                    
      }
      ktemp++; // update counter for center position

      // ============================================================
      // Loop over each of the three edges of a triangle
      // and add points along the connecting edges
      //
      // Also update the quadEdge
      // ============================================================
      ktemp3 = k3;
      ktemp4 = k4;
      ktemp5 = k5;
      ktemp6 = k6;

      for (j = 0; j < 3; j++)
      {
         // index of the middle point along any edge
         index   = g->numTriNode + g->nVertPerSide*edgeID[j]
                 + halfVertm1;

         if(strcmp(surfaceType,"sphere")==0)
            iSurface = 1;
         else 
            iSurface = 0;

         // add new vertices to the edge of a triangle      
         addVerticesOnEdge(&g->allNodePos[3*index],O,posVert,
                           &g->nodeNormal[3*index],&g->nodeNormal[3*iO],normalVert,
                           iSurface,halfVertm1,deltaDist);

         // append the position of the newly formed vertices on
         // the connecting edges to the allNodePos array         
         for (k = 0; k < halfVertm1; k++)
         {            
            g->allNodePos[3*ktemp4    ] = posVert[3*k    ];
            g->allNodePos[3*ktemp4 + 1] = posVert[3*k + 1];
            g->allNodePos[3*ktemp4 + 2] = posVert[3*k + 2];

            g->nodeNormal[3*ktemp4    ] = normalVert[3*k    ];
            g->nodeNormal[3*ktemp4 + 1] = normalVert[3*k + 1];
            g->nodeNormal[3*ktemp4 + 2] = normalVert[3*k + 2];

            tempPos[3*ktemp4  ] = g->allNodePos[3*ktemp4  ];
            tempPos[3*ktemp4+1] = g->allNodePos[3*ktemp4+1];
            tempPos[3*ktemp4+2] = g->allNodePos[3*ktemp4+2];

            ktemp4++;
         } // k loop


         if      (index == iD) // bordering Quad1 and Quad2
         {
            nodeMidDO         = ktemp4 - halfVertm1 + halfhalfVertm1; //
         }
         else if (index == iE) // bordering Quad2 and Quad3
         {
            nodeMidEO           = ktemp4 - halfVertm1 + halfhalfVertm1;
         }
         else if (index == iF) // bordering Quad1 and Quad3
         {
            nodeMidFO           = ktemp4 - halfVertm1 + halfhalfVertm1;
         }
         else
         {
            traces('Something is missing. Stopping.');
            exit(1);
         }

      } // j loop (loop over each edge)
      // ============================================================
      // Define and add the interior points and edges
      // Also define the quadConn for the interior cells using 
      // a temp array
      //
      // - Define interior points for each quad of a triangle
      // - Define interior edges of each quad of a triangle
      // - Define the quad loops for each quad of a triangle
      // ============================================================
      
      for (j = 0; j < 3; j++) // loop for each quad
      {
         if      (j == 0) // Quad ADOF
         {
            // corner points
            quadPtID[0           ][0           ] = iA;
            quadPtID[0           ][halfVertp2-1] = iF;
            quadPtID[halfVertp2-1][0           ] = iD;
            quadPtID[halfVertp2-1][halfVertp2-1] = iO;

            // other boundary points
            ktemp = halfVertm1;
            for (k = 0; k < halfVertm1; k++)
            {
               quadPtID[ktemp][0] = iD - dirAB*(k+1); // AD
               quadPtID[0][ktemp] = iF + dirCA*(k+1); // AF
               ktemp--; 

               quadPtID[halfVertp2-1][k+1] = nodeMidDO - halfhalfVertm1 + k; // DO
               quadPtID[k+1][halfVertp2-1] = nodeMidFO - halfhalfVertm1 + k; // FO
            }

         }
         else if (j == 1) // Quad DBEO
         {
            // corner points
            quadPtID[0           ][0           ] = iD;
            quadPtID[0           ][halfVertp2-1] = iO;
            quadPtID[halfVertp2-1][0           ] = iB;
            quadPtID[halfVertp2-1][halfVertp2-1] = iE;

            // other boundary points
            ktemp = halfVertm1;
            for (k = 0; k < halfVertm1; k++)
            {
               quadPtID[halfVertp2-1][ktemp] = iE - dirBC*(k+1); // BE
               quadPtID[ktemp][halfVertp2-1] = nodeMidEO - halfhalfVertm1 + k;// OE
               ktemp--; 

               quadPtID[k+1][0] = iD + dirAB*(k+1); // DB
               quadPtID[0][k+1] = nodeMidDO - halfhalfVertm1 + k;// DO
            
            }
         }
         else             // Quad OECF
         {
            // corner points
            quadPtID[0           ][0           ] = iO;
            quadPtID[0           ][halfVertp2-1] = iF;
            quadPtID[halfVertp2-1][0           ] = iE;
            quadPtID[halfVertp2-1][halfVertp2-1] = iC;

            // other boundary points
            ktemp = halfVertm1;
            for (k = 0; k < halfVertm1; k++)
            {
               quadPtID[ktemp][0] = nodeMidEO - halfhalfVertm1 + k;// OE
               quadPtID[0][ktemp] = nodeMidFO - halfhalfVertm1 + k;// OF
               ktemp--; 

               quadPtID[k+1][halfVertp2-1] = iF - dirCA*(k+1); // CF
               quadPtID[halfVertp2-1][k+1] = iE + dirBC*(k+1); // DO
            } // k loop
         } // if

         // =========================================================
         // At this point, the boundary IDs for the quads have been
         // obtained. Now obtain the interior points.
         //
         // Build middle points based on straight-lines (equal areas)
         // =========================================================
         for (jj = 0; jj < halfVertm1; jj++)
         {
            index1 = quadPtID[jj+1][0];
            index2 = quadPtID[jj+1][halfVertp2-1];
            for (kk = 0; kk < halfVertm1; kk++)
            {
               quadPtID[jj+1][kk+1] = ktemp5;

               // append into the allNodePos array
               for (k = 0; k < 3; k++)
               {
                  g->allNodePos[3*ktemp5+k] = 
                  g->allNodePos[3*index1+k] +
                  (kk+1)*deltaDist*(g->allNodePos[3*index2+k]
                                  - g->allNodePos[3*index1+k]);

                  g->nodeNormal[3*ktemp5+k] = 
                  g->nodeNormal[3*index1+k] +
                  (kk+1)*deltaDist*(g->nodeNormal[3*index2+k]
                                  - g->nodeNormal[3*index1+k]);

                  //tempPos[3*ktemp5+k] = g->allNodePos[3*ktemp5+k];

               } // k loop
               if(strcmp(surfaceType,"sphere")==0)
                  moveToBoundary(dummy1,dummy2,&g->allNodePos[3*ktemp5],
                     &g->nodeNormal[3*ktemp5],surfaceType);

               ktemp5++;
            } // kk loop
         } // jj loop

         // =========================================================
         // Build interior points based on shape-preserving ideas
         // =========================================================

         if (j == 0) // Quad ADOF
         {
            for (jj = 1; jj < halfVertp2-1; jj++)
            {
               index1 = quadPtID[jj][jj];

               // equispace the diagonals
               for (k = 0; k < 3; k++)
               {
                  g->allNodePos[3*index1+k] = g->allNodePos[3*iA+k]
                     + (jj)*deltaDist
                     *(g->allNodePos[3*iO+k] - g->allNodePos[3*iA+k]);

               } // k loop
            } // jj loop

            for (rowID = 1; rowID < halfVertp2-1; rowID++)
            {
               index2    = quadPtID[rowID][rowID];
               index3    = quadPtID[rowID][halfVertp2-1];
               deltaDistTemp = 1./(halfVertp2-1-rowID);
               counter = 1;
               for (colID = rowID+1; colID < halfVertp2-1; colID++)
               {
                  index1 = quadPtID[rowID][colID];
                  for (k = 0; k < 3; k++)
                  {
                     g->allNodePos[3*index1+k] = g->allNodePos[3*index2+k]
                        + counter*deltaDistTemp
                        *(g->allNodePos[3*index3+k] - g->allNodePos[3*index2+k]);
                  }
                  counter++;
               } // colID
            } // rowID

            for (colID = 1; colID < halfVertp2-1; colID++)
            {
               index2    = quadPtID[colID][colID];
               index3    = quadPtID[halfVertp2-1][colID];
               deltaDistTemp = 1./(halfVertp2-1-colID);
               counter = 1;

               for (rowID = colID+1; rowID < halfVertp2-1; rowID++)
               {
                  index1 = quadPtID[rowID][colID];
                   // printf("rowID,colID: [%d,%d] deltaDist: %f counter: %d\n",
                   //    rowID,colID,deltaDistTemp,counter);
                  for (k = 0; k < 3; k++)
                  {
                     g->allNodePos[3*index1+k] = g->allNodePos[3*index2+k]
                        + counter*deltaDistTemp
                        *(g->allNodePos[3*index3+k] - g->allNodePos[3*index2+k]);
                  }
                  counter++;
               } // colID
            } // rowID

         }
         else if (j == 1) // quadDOBE
         {
            for (jj = 1; jj < halfVertp2-1; jj++)
            {
               index1 = quadPtID[halfVertp2-1-jj][jj];
               // equispace the diagonals
               for (k = 0; k < 3; k++)
               {
                  g->allNodePos[3*index1+k] = g->allNodePos[3*iB+k]
                     + (jj)*deltaDist
                     *(g->allNodePos[3*iO+k] - g->allNodePos[3*iB+k]);

               } // k loop
            } // jj loop

            for (rowID = 1; rowID < halfVertp2-1; rowID++)
            {
               index2    = quadPtID[rowID][0];
               index3    = quadPtID[rowID][halfVertp2-1-rowID];
               deltaDistTemp = 1./(halfVertp2-1-rowID);
               counter = 1;
               for (colID = 1; colID < halfVertp2-1-rowID; colID++)
               {
                  // printf("rowID,colID: [%d,%d] deltaDist: %f counter: %d\n",
                  //     rowID,colID,deltaDistTemp,counter);
                  index1 = quadPtID[rowID][colID];
                  for (k = 0; k < 3; k++)
                  {
                     g->allNodePos[3*index1+k] = g->allNodePos[3*index2+k]
                        + counter*deltaDistTemp
                        *(g->allNodePos[3*index3+k] - g->allNodePos[3*index2+k]);
                  }
                  counter++;
               } // colID
            } // rowID

            for (colID = 1; colID < halfVertp2-1; colID++)
            {
               index2    = quadPtID[halfVertp2-1-colID][colID];
               index3    = quadPtID[halfVertp2-1][colID];
               deltaDistTemp = 1./colID;
               counter = 1;

               for (rowID = halfVertp2-colID; rowID < halfVertp2-1; rowID++)
               {
                  index1 = quadPtID[rowID][colID];
                   // printf("rowID,colID: [%d,%d] deltaDist: %f counter: %d\n",
                   //    rowID,colID,deltaDistTemp,counter);
                  for (k = 0; k < 3; k++)
                  {
                     g->allNodePos[3*index1+k] = g->allNodePos[3*index2+k]
                        + counter*deltaDistTemp
                        *(g->allNodePos[3*index3+k] - g->allNodePos[3*index2+k]);
                  }
                  counter++;
               } // colID
            } // rowID


         }
         else if (j == 2) // Quad OFCE
         {

            for (jj = 1; jj < halfVertp2-1; jj++)
            {
               index1 = quadPtID[halfVertp2-1-jj][halfVertp2-1-jj];

               // equispace the diagonals
               for (k = 0; k < 3; k++)
               {
                  g->allNodePos[3*index1+k] = g->allNodePos[3*iC+k]
                     + (jj)*deltaDist
                     *(g->allNodePos[3*iO+k] - g->allNodePos[3*iC+k]);

               } // k loop
            } // jj loop

            for (rowID = 1; rowID < halfVertp2-1; rowID++)
            {
               index2    = quadPtID[rowID][0];
               index3    = quadPtID[rowID][rowID];
               deltaDistTemp = 1./rowID;
               counter = 1;
               for (colID = 1; colID < rowID; colID++)
               {
                  // printf("rowID,colID: [%d,%d] deltaDist: %f counter: %d\n",
                  //     rowID,colID,deltaDistTemp,counter);
                  index1 = quadPtID[rowID][colID];
                  for (k = 0; k < 3; k++)
                  {
                     g->allNodePos[3*index1+k] = g->allNodePos[3*index2+k]
                        + counter*deltaDistTemp
                        *(g->allNodePos[3*index3+k] - g->allNodePos[3*index2+k]);
                  }
                  counter++;
               } // colID
            } // rowID


            for (colID = 1; colID < halfVertp2-1; colID++)
            {
               index2    = quadPtID[0][colID];
               index3    = quadPtID[colID][colID];
               deltaDistTemp = 1./colID;
               counter = 1;

               for (rowID = 1; rowID < colID; rowID++)
               {
                  index1 = quadPtID[rowID][colID];
                   // printf("rowID,colID: [%d,%d] deltaDist: %f counter: %d\n",
                   //    rowID,colID,deltaDistTemp,counter);
                  for (k = 0; k < 3; k++)
                  {
                     g->allNodePos[3*index1+k] = g->allNodePos[3*index2+k]
                        + counter*deltaDistTemp
                        *(g->allNodePos[3*index3+k] - g->allNodePos[3*index2+k]);
                  }
                  counter++;
               } // colID
            } // rowID
         }

         // =========================================================
         // At this point, the boundary IDs for the quads have been
         // obtained. Now obtain the interior points.
         //
         // Build middle points based on straight-lines (equal areas)
         // =========================================================
         
         for (jj = 0; jj < halfVertm1; jj++)
         {
            index1 = quadPtID[jj+1][0];
            index2 = quadPtID[jj+1][halfVertp2-1];
            for (kk = 0; kk < halfVertm1; kk++)
            {
               quadPtID[jj+1][kk+1] = ktemp6;

               // append into the allNodePos array
               for (k = 0; k < 3; k++)
               {
                  tempPos[3*ktemp6+k] = 
                  g->allNodePos[3*index1+k] +
                  (kk+1)*deltaDist*(g->allNodePos[3*index2+k]
                                  - g->allNodePos[3*index1+k]);

                  g->nodeNormal[3*ktemp6+k] = 
                  g->nodeNormal[3*index1+k] +
                  (kk+1)*deltaDist*(g->nodeNormal[3*index2+k]
                                  - g->nodeNormal[3*index1+k]);

                  //tempPos[3*ktemp5+k] = g->allNodePos[3*ktemp5+k];

               } // k loop
               if(strcmp(surfaceType,"sphere")==0)
                  moveToBoundary(dummy1,dummy2,&g->allNodePos[3*ktemp6],
                     &g->nodeNormal[3*ktemp6],surfaceType);

               ktemp6++;
            } // kk loop
         } // jj loop


      } // j loop (loop over each quad)
      k1++;         // for center
      k2 += 3*pow4; // for number of cells in a triangle
      k3  = ktemp3; // edges inside triangle
      k4  = ktemp4; // vertices on int edges
      k5  = ktemp5; // quad int points
      k6  = ktemp6; // quad int points

	} // i loop (g->numTriangle)


   // recombine based on weights (1 = straight, 0 = shape-preserving)
   weight = 0.5;
   for (i = 0; i < threenum; i++)
      g->allNodePos[i] = weight*tempPos[i] + (1.-weight)*g->allNodePos[i];


   free(tempPos);
}

// ##################################################################
// END OF FILE
// ##################################################################
