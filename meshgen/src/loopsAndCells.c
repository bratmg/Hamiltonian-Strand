// ##################################################################
//
// loopsAndCells.c
//
// Routines related to the generation of 
// triangular/quad loops and quad cells
// ##################################################################
#include <stdio.h>
#include <stdlib.h>

#include "globalVariables.h"
#include "meshtype.h"
#include "loopsAndCells.h"

// ##################################################################
// 
// Find triangular and vertex loops
//
// ##################################################################
void triCellAndVertexLoops(GRID *g)
{

   printf("#meshgen: Creating triangle and vertex loops ...\n");

   int i,i1,i2,j,j1,k,kk,k1,k2;
   int ikc,ikv,icont;
   int ncl,nc,nv;
   int V1,V2;
   int twoj,twojp1,threej,threejp1,threejp2;
   int tempLen,nodeID;
   int nact,kloop1,kloop2;

   // Allocation for the loop and initialization of the same
   g->vertLoop = (LOOP *) malloc(sizeof(LOOP));
   g->triLoop  = (LOOP *) malloc(sizeof(LOOP));

	// identify max length of a loop (useful for mem allocation)
   g->triLoop->totLen   = g->tri2nodeIndex[g->numTriNode];
   g->vertLoop->totLen  = g->triLoop->totLen + g->numTriNode;

   g->triLoop->index    = (int *) malloc(sizeof(int)*(g->numTriNode + 1));
   g->vertLoop->index   = (int *) malloc(sizeof(int)*(g->numTriNode + 1));

   g->triLoop->ID       = (int *) malloc(sizeof(int)*g->triLoop->totLen);
   g->vertLoop->ID      = (int *) malloc(sizeof(int)*g->vertLoop->totLen);

	// Maximum sizes of triangular and vertex loops
   g->triLoop->maxLen= -1;
   for (i = 1; i < g->numTriNode + 1; i++)
   {
      tempLen = g->tri2nodeIndex[i] - g->tri2nodeIndex[i-1];
      
      if (tempLen > g->triLoop->maxLen)
         g->triLoop->maxLen = tempLen;
   }

   g->vertLoop->maxLen = g->triLoop->maxLen+1;

   // initialize the tri and vert loop indices to 0
   for (i = 0; i < g->numTriNode + 1; i++)
      g->triLoop->index[i] = 0;

   for (i = 0; i < g->numTriNode + 1; i++)
      g->vertLoop->index[i] = 0;

   // Initialize cell connectivity array (size = maxLen*3 and maxLen*2)
   int CC[3*g->triLoop->maxLen];
   int CC2[2*g->triLoop->maxLen];
   int cLoop1[g->triLoop->maxLen];
   int cLoop2[g->triLoop->maxLen];
   int iact[g->triLoop->maxLen]; 
   int EV[2*g->vertLoop->maxLen];

   // Loop through the total number of nodes
   ikc = 0; 
   ikv = 0;
	for (i = 0; i < g->numTriNode; i++)
	{

		k  = g->triLoop->maxLen;
      k1 = k;
      k2 = k+3;

      // number of triangles around a given node
      i1 = g->tri2nodeIndex[i];     // start index for a given node
      i2 = g->tri2nodeIndex[i+1]-1; // end index for a given node
      nc = i2-i1+1;                 // number of cells

		// Creates the list of cell connectivities
		// e.g.: triangles connected to node 1
 	   //        128 100  1
 	   //        157 128  1
 	   //        157   1  2
      j1 = 0;
      for (j = i1; j <= i2; j++ )
      {
         CC[j1]   = g->triConn[3*(g->tri2nodeList[j])];
         CC[j1+1] = g->triConn[3*(g->tri2nodeList[j]) + 1];
         CC[j1+2] = g->triConn[3*(g->tri2nodeList[j]) + 2];         

         j1 += 3;
      } // j loop
            
      // create the CC2 array
      // Contains the OTHER two node IDs for a given triangle
      // Loop for each of the connected cells of a node   
      for (j = 0; j < nc ; j++)
      {
         twoj     = 2*j;
         twojp1   = twoj+1;
         threej   = 3*j;
         threejp1 = threej+1;
         threejp2 = threej+2;

         
         if( CC[threej] == i)         
         {
            // If the 1st point of three is the vertex itself
            // the other two are 2 and 3
            CC2[twoj]   = CC[threejp1];
            CC2[twojp1] = CC[threejp2];
         }
         else if (CC[threejp1] == i)
         {
            // If the 2nd point of three is the vertex itself
            // the other two are 3 and 1
            CC2[twoj]   = CC[threejp2];
            CC2[twojp1] = CC[threej];
         }
         else
         {
            // If the 3rd point of three is the vertex itself
            // the other two are 1 and 2
            CC2[twoj]   = CC[threej];
            CC2[twojp1] = CC[threejp1];
         }


      } // j loop
      

      // If the number of connections is only 1
      // (Haven't encountered this condition yet)
      if (nc == 1)
      {
         ncl = 1; // number of cells
         nv  = 2; // number of vertices
         g->triLoop->ID[ikc]    = g->tri2nodeList[i1];
         g->vertLoop->ID[ikv]   = CC2[0];
         g->vertLoop->ID[ikv+1] = CC2[1];

         // update counters
         ikc += ncl;
         ikv += nv;

      }
      else // if more than 1 triangle to a node (more common)
      {

         // cLoop = cell loop
         cLoop2[0] = g->tri2nodeList[i1]; // first element of the node itself
         kloop1    = -1;                   // counter for number of elements in cLoop1
         kloop2    =  0;                   // counter for number of elements in cLoop2

         V1        = CC2[0];
         V2        = CC2[1];
         EV[k+1]     = V1; // end vertices for the loop
         EV[k+2]   = V2; // the last two vertices are CC2[0] and CC2[1]

         // set iact to 1s
         for (kk = 0; kk < g->triLoop->maxLen; kk++)
            iact[kk] = 1;

         iact[0]   = 0; 
         nact      = nc-1; // number of active edges connections
         icont     = 1;    // flag for continuation

         // while continuation flag is 1
         // This algorithm builds the vertices from a list of
         // CC values. The chain builds to left if V1 and to the right if V2
         // For e.g., stop at i == 906 (look for ppt in the folder explaining
         // this algorithm)
         while (icont == 1)
         {
            // loop from 1 to the remaining number cells
            for (j = 1; j < nc; j++)
            {
               // if the connection is still active
               if (iact[j] == 1)
               {
                  // loop for the other two vertices
                  for (kk = 0; kk < 2; kk++)
                  {
                     // The (1-kk) basically switches to the other column
                     // of the same row in CC2. If kk = 1, switches to 0,
                     // is kk = 0, switches to 1
                     if (CC2[2*j+kk] == V1)
                     {                      
                        V1               = CC2[2*j+1-kk]; // switch V1 other column in same row
                        EV[k1]           = V1;            // fill EV from element maxLen to 1
                        iact[j]          = 0;
                        cLoop1[kloop1+1] = g->tri2nodeList[i1+j]; // append nodeID
                        k1--;                        
                        nact--;                        
                        kloop1++;
                        break;
                     }
                     else if (CC2[2*j+kk] == V2)
                     {                        
                        V2               = CC2[2*j+1-kk];
                        EV[k2]           = V2;
                        iact[j]          = 0;
                        cLoop2[kloop2+1] = g->tri2nodeList[i1+j]; // append nodeID
                        k2++;                        
                        nact--;                        
                        kloop2++;
                        break;
                     }


                  } // kk for loop                  
               } // if iact


               // Termination condition for icont based on if there are
               // any remaining active edges left (nact)
               if (nact == 0)
               {
                  icont = 0;
                  break;
               }

            } // j for loop            
         } // while condition (icont)
      
         // cell loops
         ncl = kloop1 + kloop2 + 2; // number of cells         

         for (kk = kloop1; kk >= 0; kk--)
         {
            g->triLoop->ID[ikc] = cLoop1[kk];
            ikc++;
         }

         for (kk = 0; kk < kloop2+1; kk++)
         {            
            g->triLoop->ID[ikc] = cLoop2[kk];
            ikc++;
         }

         // vertex loops
         nv = k2-k1-1; // number of vertices         
         for (kk = k1+1; kk < k2; kk++)
         {
            g->vertLoop->ID[ikv] = EV[kk];
            ikv++;
         }
      } // if nc == 1

      // create the indices of the loops
      g->triLoop->index[i+1]  = g->triLoop->index[i]  + ncl;
      g->vertLoop->index[i+1] = g->vertLoop->index[i] + nv;

	} // i loop
   
}

// ##################################################################
// 
// Find inner and outer quad loops
//
// ##################################################################
void quadLoops(GRID *g)
{
   
   printf("#meshgen: Creating quadrilateral loops ...\n");

   int i,j,k,ii,jj,kk,kk1,kk2,kk3;
   int iA,iB,iC,ncl,is2,tempIdx,offset;
   int edgeID[3],edge[2];
   int cell,cell3;
   int *CL, *qeVec1, *qeVec2, *qeVec3, **qeVec;
   int *q1loop, *q2loop, *iqloop, **qloop;

   // metrics used for edges and vertices
   int numEdgep1       = g->nEdgePerSide+1;
   int halfEdgePerSide = 0.5*g->nEdgePerSide;
   int halfVertm1      = 0.5*(g->nVertPerSide-1);
   int halfhalfVertm1  = 0.5*(halfVertm1-1);
   int halfVertp2      = halfVertm1 + 2;
   int horzEdgePerQuad = halfEdgePerSide*halfVertm1;
   int vertEdgePerQuad = halfEdgePerSide*halfVertm1;
   int edgePerQuad     = horzEdgePerQuad+vertEdgePerQuad;
   int edgePerLine     = halfEdgePerSide;

   iqloop = (int *) malloc(sizeof(int)*(g->numTriNode+1));
   
   // allocate memory
   CL     = (int *) malloc(sizeof(int)*g->triLoop->maxLen);
   iqloop = (int *) malloc(sizeof(int)*(g->numTriNode+1));
   
   // allocate memory for the loops
   tempIdx = numEdgep1*g->triLoop->maxLen+1;
   qeVec  = (int **) malloc(sizeof(int *)*pow2);
   for (i = 0; i < pow2; i++)
      qeVec[i]  = (int *) malloc(sizeof(int)*tempIdx);

   for (i = 0; i < pow2; i++)
      for (j = 0; j < tempIdx; j++)
         qeVec[i][j] = -1;

   tempIdx = numEdgep1*g->numTriNode*g->triLoop->maxLen + 1;
   qloop  = (int **) malloc(sizeof(int *)*pow2);
   for (i = 0; i < pow2; i++)
      qloop[i] = (int *) malloc(sizeof(int)*tempIdx);

   for (i = 0; i < pow2; i++)
      for (j = 0; j < tempIdx; j++)
         qloop[i][j] = -1;


   for (i = 0; i < g->numTriNode+1; i++ )
      iqloop[i] =  0;

   // Running counters for quad loops
   kk = 0; kk1 = 0; kk2 = 0; kk3 = 0;

   //
   // Loop through the total number of nodes
   for (ii = 0; ii < g->numTriNode; ii++)
   {
      i   = g->colourIndex[ii];         

      k = 0; // counter for number of cells in loop
      for (j = g->triLoop->index[i]; j < g->triLoop->index[i+1]; j++)
      {
         CL[k] = g->triLoop->ID[j];      
         k++;         
      }
      
      ncl = k; // number of cells in loop

      // Loop through the triangular cells in the loop
      // Refer documentation
      for (j = 0; j < ncl; j++)
      {
      
         cell  = CL[j];
         cell3 = 3*cell;        

         iA = g->triConn[cell3  ];
         iB = g->triConn[cell3+1];
         iC = g->triConn[cell3+2];      

         // list of the edge ID for a given triangle
         edgeID[0] = g->edge2triList[cell3  ];
         edgeID[1] = g->edge2triList[cell3+1];
         edgeID[2] = g->edge2triList[cell3+2];      
      
         // =========================================================
         // if node A ID is the colour index of cell (group it)
         // =========================================================
         if (iA == i)
         {
            // traces( ---------- NODE A -----------);
            // traces(maybe 3*.. to cell*...);exit(1);
            // Loop through the edges of the triangle
            for (k = 0; k < 3; k++)
            {
               // list columns 1--4 of a particular edge
               edge[0] = g->triEdge[edgeID[k]][0];
               edge[1] = g->triEdge[edgeID[k]][1];               

               is2     = g->nEdgePerSide*edgeID[k]; // index of edge in Qedge array

               // edge A-B
               if (edge[0] == iA && edge[1] == iB)
               {
                  // traces(AB);
                  for (jj = 0; jj < pow2; jj++) // number of loops
                     qeVec[jj][numEdgep1*j] = is2 + halfEdgePerSide + jj;
               }

               // edge B-A
               if (edge[0] == iB && edge[1] == iA)
               {
                  // traces(BA);
                  for (jj = 0; jj < pow2; jj++)
                     qeVec[jj][numEdgep1*j] = is2 + halfEdgePerSide - (jj+1);
               }

               // edge C-A
               if (edge[0] == iC && edge[1] == iA)
               {
                  // traces(CA);
                  for (jj = 0; jj < pow2; jj++)
                     qeVec[jj][numEdgep1*j + g->nEdgePerSide] = is2 + halfEdgePerSide - (jj+1);
               }

               // edge A-C
               if (edge[0] == iA && edge[1] == iC)
               {
                  // traces(AC);
                  for (jj = 0; jj < pow2; jj++)
                     qeVec[jj][numEdgep1*j + g->nEdgePerSide] = is2 + halfEdgePerSide + jj;
               }


               // collect vertical edges in quads 2 and 3
               for (jj = 0; jj < pow2; jj++) // number of loops
               {
                  // ================================================
                  // QUAD 2
                  // ================================================
                  for (kk = 0; kk < halfVertm1; kk++) // along each loop
                  {
                     offset  = g->numTriEdge*g->nEdgePerSide 
                             + 3*halfEdgePerSide*(cell+1)
                             + cell*3*edgePerQuad
                             + edgePerQuad + horzEdgePerQuad;

                     tempIdx =  offset + jj + kk*edgePerLine;

                     qeVec[jj][numEdgep1*j + (kk+1)] = tempIdx;
                  }
                  // ================================================
                  // LINE EO
                  // ================================================
                  offset  = g->midEdgeID[3*cell+1] + halfEdgePerSide;
                  
                  tempIdx = offset - (jj+1);

                  qeVec[jj][numEdgep1*j + halfVertp2-1] = tempIdx;
                  // 
                  // QUAD 3
                  //
                  for (kk = 0; kk < halfVertm1; kk++) // along each loop
                  {
                     offset  = g->numTriEdge*g->nEdgePerSide 
                             + 3*halfEdgePerSide*(cell+1) 
                             + cell*3*edgePerQuad
                             + 2*edgePerQuad + horzEdgePerQuad;

                     tempIdx =  offset + jj + kk*edgePerLine;

                     qeVec[jj][numEdgep1*j + (kk+halfVertp2)] = tempIdx;
                  }

                  // if (1)
                  // {
                  //    for (kk = 0; kk < numEdgep1; kk++)
                  //       printf("(%d)(%d) [] %d\n",jj,kk,qeVec[jj][kk]);
                  // }

               } // jj loop
               // traces(END iA================);               


            } // k loop

         } // if 
         // =========================================================
         // If NODE B is of the same colour
         // =========================================================
         else if ( iB == i)
         {
            // traces( ---------- NODE B -----------);
            // Loop through the edges of the triangle
            for (k = 0; k < 3; k++)
            {
               // list columns 1--4 of a particular edge
               edge[0] = g->triEdge[edgeID[k]][0];
               edge[1] = g->triEdge[edgeID[k]][1];  

               is2     = g->nEdgePerSide*edgeID[k]; // index of edge in Qedge array


               // edge B-C
               if (edge[0] == iB && edge[1] == iC)
               {
                  // traces(BC);
                  for (jj = 0; jj < pow2; jj++) // number of loops
                     qeVec[jj][numEdgep1*j] = is2 + halfEdgePerSide + jj;
               }

               // edge C-B
               if (edge[0] == iC && edge[1] == iB)
               {
                  // traces(CB);
                  for (jj = 0; jj < pow2; jj++)
                     qeVec[jj][numEdgep1*j] = is2 + halfEdgePerSide - (jj+1);
               }

               // edge A-B
               if (edge[0] == iA && edge[1] == iB)
               {
                  // traces(AB);
                  for (jj = 0; jj < pow2; jj++)
                     qeVec[jj][numEdgep1*j + g->nEdgePerSide] = is2 + halfEdgePerSide - (jj+1);
               }

               // edge B-A
               if (edge[0] == iB && edge[1] == iA)
               {
                  // traces(BA);
                  for (jj = 0; jj < pow2; jj++)
                     qeVec[jj][numEdgep1*j + g->nEdgePerSide] = is2 + halfEdgePerSide + jj;
               }

               // collect vertical edges in quads 3 and 1
               for (jj = 0; jj < pow2; jj++) // number of loops
               {
                  // ================================================
                  // QUAD 3
                  // ================================================
                  for (kk = 0; kk < halfVertm1; kk++) // along each loop
                  {
                     offset  = g->numTriEdge*g->nEdgePerSide 
                             + 3*halfEdgePerSide*(cell+1) 
                             + cell*3*edgePerQuad
                             + 2*edgePerQuad;

                     tempIdx =  offset + (halfVertm1-1)*halfEdgePerSide + jj
                             - kk*halfEdgePerSide;

                     qeVec[jj][numEdgep1*j + (kk+1)] = tempIdx;
                  }
                  // ================================================
                  // LINE FO
                  // ================================================
                  offset  = g->midEdgeID[3*cell+2] + halfEdgePerSide;
                  
                  tempIdx = offset - (jj+1);

                  qeVec[jj][numEdgep1*j + halfVertp2-1] = tempIdx;
                  // ================================================
                  // QUAD 1
                  // ================================================
                  for (kk = 0; kk < halfVertm1; kk++) // along each loop
                  {
                     offset  = g->numTriEdge*g->nEdgePerSide 
                             + 3*halfEdgePerSide*(cell+1) 
                             + cell*3*edgePerQuad
                             + 0*edgePerQuad + horzEdgePerQuad;

                     tempIdx =  offset + vertEdgePerQuad 
                             - (jj+1) - kk*halfEdgePerSide;

                     qeVec[jj][numEdgep1*j + (kk+halfVertp2)] = tempIdx;
                  }


               }
               // traces(END iB================);
            } // k loop
            

          } // if condition
          // ========================================================
          // NODE C
          // ========================================================
          else if (iC == i)
          {
             // traces( ---------- NODE C -----------);
            // Loop through the edges of the triangle
            for (k = 0; k < 3; k++)
            {
               // list columns 1--4 of a particular edge
               edge[0] = g->triEdge[edgeID[k]][0];
               edge[1] = g->triEdge[edgeID[k]][1];               

               is2     = g->nEdgePerSide*edgeID[k]; // index of edge in Qedge array

               // edge C-A
               if (edge[0] == iC && edge[1] == iA)
               {
                  // traces(CA);trace(is2);trace(halfEdgePerSide);
                  for (jj = 0; jj < pow2; jj++) // number of loops
                     qeVec[jj][numEdgep1*j] = is2 + halfEdgePerSide + jj;
               }

               // edge A-C
               if (edge[0] == iA && edge[1] == iC)
               {
                  // traces(AC);
                  for (jj = 0; jj < pow2; jj++)
                     qeVec[jj][numEdgep1*j] = is2 + halfEdgePerSide - (jj+1);
               }

               // edge B-C
               if (edge[0] == iB && edge[1] == iC)
               {
                  // traces(BC);
                  for (jj = 0; jj < pow2; jj++)
                     qeVec[jj][numEdgep1*j + g->nEdgePerSide] = is2 + halfEdgePerSide - (jj+1);
               }

               // edge C-B
               if (edge[0] == iC && edge[1] == iB)
               {
                  // traces(AC);
                  for (jj = 0; jj < pow2; jj++)
                     qeVec[jj][numEdgep1*j + g->nEdgePerSide] = is2 + halfEdgePerSide + jj;
               }

               // collect vertical edges in quads 3 and 1
               for (jj = 0; jj < pow2; jj++) // number of loops
               {
                  // ================================================
                  // QUAD 1
                  // ================================================
                  for (kk = 0; kk < halfVertm1; kk++) // along each loop
                  {
                     offset  = g->numTriEdge*g->nEdgePerSide 
                             + 3*halfEdgePerSide*(cell+1) 
                             + cell*3*edgePerQuad
                             + 0*edgePerQuad;

                     tempIdx =  offset + halfEdgePerSide - (jj+1)
                             + kk*halfEdgePerSide;

                     qeVec[jj][numEdgep1*j + (kk+1)] = tempIdx;
                  }
                  // ================================================
                  // LINE DO
                  // ================================================
                  offset  = g->midEdgeID[3*cell] + halfEdgePerSide;
                  tempIdx = offset - (jj+1);

                  qeVec[jj][numEdgep1*j + halfVertp2-1] = tempIdx;
                  // ================================================
                  // QUAD 2
                  // ================================================
                  for (kk = 0; kk < halfVertm1; kk++) // along each loop
                  {
                     offset  = g->numTriEdge*g->nEdgePerSide 
                             + 3*halfEdgePerSide*(cell+1) 
                             + cell*3*edgePerQuad
                             + 1*edgePerQuad;

                     tempIdx =  offset + halfEdgePerSide - (jj+1)
                             + kk*halfEdgePerSide;

                     qeVec[jj][numEdgep1*j + (kk+halfVertp2)] = tempIdx;
                  }

                  // for (kk = 0; kk < numEdgep1; kk++)
                  //    trace(qeVec[jj][kk]);

               }
               // traces(END iC================);

            } // k loop

          }
          else
          {            
            traces(Something is erroneous. Stopping.);
            exit(1);
         }

      } // j loop
      
      // ============================================================
      // Build full quad loops (without repitition)
      // qeVec are local arrays for a given node
      // qloops are accumulative loops for all nodes,
      // which are indexed based on quadLoop->index. Single index
      // for both are sifficient as Q1 and Q2 are of the same length
      // ============================================================

      // ============================================================
      // qloop for the loops
      // ============================================================
      for (jj = 0; jj < pow2; jj++) // number of loops
         qloop[jj][kk1] = qeVec[jj][0];

      for (j = 0; j < ncl; j++)
      {
         for (jj = 0; jj < pow2; jj++)
         {
            for (k = 1; k < numEdgep1; k++)
            {
               qloop[jj][kk1 + k] = qeVec[jj][numEdgep1*j + k];
               // printf("(%d)(%d) [qloop] %d\n",jj,kk1+k,qloop[jj][kk1+k]);
            }
         }
         kk1 += g->nEdgePerSide;

      } // jj (number of loops) 
      kk1++;
      //traces('stopping here');exit(1);

      // start and end index of qloop for a given node
      iqloop[ii+1] = iqloop[ii] + g->nEdgePerSide*ncl + 1;
   } // ii loop
   
   // ===============================================================
   // Copy the loops into the quadLoop struct
   // ===============================================================
   int i1,i2,nel;
   int nloops = pow2*kk1;

   // allocate memory for quadLoop and initialization
   g->quadLoop = (LOOP *) malloc(sizeof(LOOP));

   g->quadLoop->ID    = (int *) malloc(sizeof(int)*nloops);
   for (i = 0; i < nloops; i++)
      g->quadLoop->ID[i] = 0;

   g->quadLoop->index = (int *) malloc(sizeof(int)*(pow2*g->numTriNode+1));
   for (i = 0; i < 2*g->numTriNode+1 ; i++)
      g->quadLoop->index[i] = 0;
   
   g->quadLoop->maxLen = -1;

   // reset counters
   k = 0; kk = 0;

   // Run through all the nodes
   // quadLoop->ID are ordered as
   // [q1loop (node1) | q2loop (node1) | q1loop (node2) | q2loop (node2) ... ]
   for (i = 0; i < g->numTriNode; i++)
   {
      i1  = iqloop[i];   // start index for the loop
      i2  = iqloop[i+1]; // end index for the loop
      nel = i2-i1;       // total number of edges for the loop
      // set the maximum length of a quad loop
      g->quadLoop->maxLen = MAX(g->quadLoop->maxLen,nel);

      // ordering of loop IDs
      for (jj = 0; jj < pow2; jj++) 
      {
         for (j = 0; j < nel; j++)
         {
            g->quadLoop->ID[k+j] = qloop[jj][i1+j];
             // printf("(%d)(%d) [qloop] %d\n",jj,i1+j,qloop[jj][i1+j]);
         }
         k += nel;

         // ordering of loop indices
         g->quadLoop->index[kk+1] = g->quadLoop->index[kk] + nel;
         kk++;
      }
   }

   // set total length of the loops
   g->quadLoop->totLen = k;

}

// ##################################################################
// END OF FILE
// ##################################################################
