// ##################################################################
//
// edgeOperations.c
// 
// Associated edge operations
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "globalVariables.h"
#include "meshtype.h"
#include "flushToSurface.h"
#include "smoothOperations.h"
#include "edgeOperations.h"

#define N_EDGE 4
#define N_QEDGE 6
#define ONE_THIRD 0.3333333333333333


#define INDX( row, col, ld) ( row*ld + col )
// ##################################################################
//
// findTriangleEdges
//
// Function to find the edges of a polygon. Information of an edge
// is outlined as the two node IDs of an edge, the triangle ID it
// belongs to and the triangle ID the edge is shared with (if any)
// ##################################################################
void findTriangleEdges(GRID *g)
{

   printf("#meshgen: Finding edges on the triangle ...\n");

   int   i,j;
   int   me,ne;
   int  *iPtr,ip; 
   int **eTmp,eLoc[2],e1[2],e2[2];
   int   iFlag;

   int   halfEdgePerSide = 0.5*g->nEdgePerSide;
   int   halfVertm1      = 0.5*(g->nVertPerSide-1);

   me = g->nEdgePerSide*g->numTriangle; // maximum number of edges (just a safe upper bound?)

   iPtr  = (int *) malloc(sizeof(int)*g->numTriNode);
   for (i = 0; i < g->numTriNode; i++)
      iPtr[i] = -1;

   // allocate of **eTmp
   eTmp    = (int **) malloc(sizeof(int *)*me);
   for (i = 0; i < me; i++)
      eTmp[i] = (int *) malloc(sizeof(int)*(N_EDGE+1));

   // Initialize to -1
   for (i = 0; i < me; i++)
      for (j = 0; j <+ N_EDGE; j++)
         eTmp[i][j] = -1;

   ne = 0; // set in insertEdge function

   // loop over total number of cells
   for (i = 0; i < g->numTriangle; i++)
   {
      // loop over number of vertices per triangle
      for (j = 0; j < 3; j++)
      {
         
         iFlag   = 1;

         eLoc[0] = g->triConn[ 3*i +  j       ];
         eLoc[1] = g->triConn[ 3*i + ((j+1) % 3)];            

         // create temporary variable
         e1[0]   = eLoc[0];
         e1[1]   = eLoc[1];

         // swap procedure to ensure e1[1] > e1[0]
         if (e1[0] > e1[1])
            SWAP(e1[0],e1[1]);

         ip = iPtr[e1[0]];         

         // section to identify common edges
         while (ip > -1 && iFlag == 1)
         {
            e2[0] = eTmp[ip][0];
            e2[1] = eTmp[ip][1];

            // swap procedure to ensure e2[1] > e2[0]
            if (e2[0] > e2[1])
               SWAP(e2[0],e2[1]);

            // if both edges e1 and e2 are identical
            if (abs(e1[0]-e2[0]) == 0 && abs(e1[1]-e2[1]) == 0 )
            {               
               eTmp[ip][3] = i;
               iFlag       = 0;
            }

            ip = eTmp[ip][4];

         } // while loop

         if (iFlag == 1)
         {            
            eTmp[ne][0] = eLoc[0];
            eTmp[ne][1] = eLoc[1];
            eTmp[ne][2] = i;
            eTmp[ne][4] = iPtr[e1[0]];
            iPtr[e1[0]] = ne;
            ne++;
         }


      } // j loop
   } // i loop   


	// allocate g->triEdge->ID (hopefully ne is less than me)
   g->triEdge = (int **) malloc(sizeof(int *)*ne);
   for (i = 0; i < ne; i++)
   	g->triEdge[i] = (int *) malloc(sizeof(int)*N_EDGE);

	// 4 element block
	// 1st and 2nd elements are the node IDs of the edge
	// 3rd element is the trinagle ID the edge belongs to
	// 4th element is the triangle ID with which the edge is shared (0 is none)
	for (i = 0; i < ne; i++)
	{		
		g->triEdge[i][0] = eTmp[i][0];
		g->triEdge[i][1] = eTmp[i][1];
		g->triEdge[i][2] = eTmp[i][2];
		g->triEdge[i][3] = eTmp[i][3];

      if (nGrid > 1)
      {
         for (j = 0; j < g->boundaryCount; j++)
            if (i == g->boundaryID[j])
               g->triEdge[i][3] = eTmp[i][3] = -5;
      }
	} // i loop

	g->numTriEdge = ne; 


   // nodePts arrays contains the list of all points in the domain
   // number of points =   number of triangle nodes + (original points)
   //                    3*number of triangle edges + (those divided at edge)
   //                    7*number of interior points  (of the triangle)
   // Although separate arrays are created for each, it is important to 
   // consolidate all node points into one as the loop indexing is based on
   // the single array   
   g->numNodePos = g->numTriNode 
                 + g->nVertPerSide*g->numTriEdge 
                 + g->numTriangle*(3*halfVertm1 + 1 + 3*halfVertm1*halfVertm1);
   g->allNodePos  = (double *) malloc(sizeof(double)*3*g->numNodePos);
   g->surfNodePos = (double *) malloc(sizeof(double)*3*g->numNodePos);

   // append original node points
   j = 0;
   for (i = 0; i < g->numTriNode; i++)
   {
      g->allNodePos[3*i  ] = g->nodePosTri[j  ];
      g->allNodePos[3*i+1] = g->nodePosTri[j+1];
      g->allNodePos[3*i+2] = g->nodePosTri[j+2];
      j += 3;
   }

   // Free memory used
   free(iPtr);
   
   //for (i = 0; i < me; i++);
   //   free(eTmp[i]);   
   free(eTmp);

}

// ##################################################################
//
// createVerticesOnEdge
//
// Function to create vertices on the triangle edges to form the 
// basis for the different quad cells
//
// ##################################################################
void createVerticesOnEdge(GRID *g)
{

   printf("#meshgen: Creating vertices on triangular edges ...\n");

   int     i,i1,i2,j,k1,k2,ktemp;
   int     iSurface,nVert,nEdge;
   int    *vertID;
   double  posA[3],posB[3],normalA[3],normalB[3],deltaDist;
   double *posVert,*normalVert;
   int     halfEdgePerSide = 0.5*g->nEdgePerSide;
   int     halfVertm1      = 0.5*(g->nVertPerSide-1);

   k1 = g->numTriNode;
   k2 = 0;

   // allocate
   int n = g->nEdgePerSide*g->numTriEdge// + 25*g->numTriangle;
         + g->numTriangle*(3*halfEdgePerSide + 3*(2*halfVertm1*halfEdgePerSide));
   g->quadEdge = (int **) malloc(sizeof(int *)*n);
   for (i = 0; i < n; i++)
      g->quadEdge[i] = (int *) malloc(sizeof(int)*(N_QEDGE));

   for (i = 0; i < n; i++)
      for (j = 0; j < N_QEDGE; j++)
         g->quadEdge[i][j] = -1;

   // create unique cells ID for edges that border a subdomain
   if (nGrid > 1) subdomainEdgeBlanking(g);

   g->iBoundary = (int *) malloc(sizeof(int)*g->numNodePos);
   for (i = 0; i < g->numNodePos; i++)
      g->iBoundary[i] = 0;

   // allocate posVert   
   deltaDist  = (double) 1./g->nEdgePerSide;
   posVert    = (double *) malloc(sizeof(double)*3*g->nVertPerSide);
   normalVert = (double *) malloc(sizeof(double)*3*g->nVertPerSide);

   vertID    = (int *) malloc(sizeof(int)*(g->nVertPerSide+2));

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
         posA[j]    = g->allNodePos[3*i1 + j];
         posB[j]    = g->allNodePos[3*i2 + j];  
         normalA[j] = g->nodeNormal[3*i1 + j];
         normalB[j] = g->nodeNormal[3*i2 + j];
      }

      // Is this a domain specific condition. Must check.
      if ( strcmp(surfaceType,"naca")==0 
            && (g->triEdge[i][3]==-1) 
            && (abs(posA[1])<0.61) && (abs(posA[0])<1.1) )
         iSurface = 1;                  
      else if (strcmp(surfaceType,"sphere")==0 ||
               strcmp(surfaceType,"robin")==0 )
         iSurface = 1;
      else
         iSurface = 0;      

      // add new vertices to the edge of a triangle      
      addVerticesOnEdge(posA,posB,posVert,
                        normalA,normalB,normalVert,iSurface,
                        g->nVertPerSide,deltaDist);

      if(iSurface==1)
      {
         g->iBoundary[i1  ] = 1; // posA
         g->iBoundary[i2  ] = 1; // posB
         g->iBoundary[k1  ] = 1; // posA1
         g->iBoundary[k1+1] = 1; // posA2
         g->iBoundary[k1+2] = 1; // posA3
      }

      // append the position of the newly formed vertices onto
      // the edge array, and update their normals as well
      ktemp = k1;       
      for (j = 0; j < g->nVertPerSide; j++)
      {
         g->allNodePos[3*ktemp    ] = posVert[3*j    ];
         g->allNodePos[3*ktemp + 1] = posVert[3*j + 1];
         g->allNodePos[3*ktemp + 2] = posVert[3*j + 2];

         g->nodeNormal[3*ktemp    ] = normalVert[3*j    ];
         g->nodeNormal[3*ktemp + 1] = normalVert[3*j + 1];
         g->nodeNormal[3*ktemp + 2] = normalVert[3*j + 2];

         ktemp++;
      } // j loop

      // specify quadedges
      //vertID = [i1 k1 k1+1 k1+2 ... k1+(nVertPerSide-1) i2]
      vertID[0] = i1; vertID[g->nVertPerSide+1] = i2;

      for (j = 0; j < g->nVertPerSide; j++)
         vertID[j+1] = k1+j;

      // The following incides are a one-to-one map to the indices
      // of the edges and nodes to the g->allNodePos
      // Qedges(k2+1) [ A1 - A  ]
      // Qedges(k2+2) [ A2 - A1 ]
      // Qedges(k2+3) [ A2 - A3 ]
      // Qedges(k2+4) [ A3 - B  ]
      ktemp = k2; 
      for (j = 0; j < g->nEdgePerSide; j++)
      {
         if (j < 0.5*g->nEdgePerSide)
         {
            g->quadEdge[ktemp][0] = vertID[j+1];
            g->quadEdge[ktemp][1] = vertID[j];
         }
         else
         {
            g->quadEdge[ktemp][0] = vertID[j];
            g->quadEdge[ktemp][1] = vertID[j+1];
         }
         ktemp++;

      } // j loop

      k1 += g->nVertPerSide; 
      k2 += g->nEdgePerSide; 

   } // i loop

   // total number fo elements in Q 
   // (hopefully lesser than 4*nEt + 18*nt)
   g->numQuadEdge = k2;

   free(posVert);


   // for (i = 0 ; i < k2; i++)
      // printf("quadEdge[i][0]->quadEdge[i][1]: %d -> %d\n",g->quadEdge[i][0],g->quadEdge[i][1]);


}

// ##################################################################
//
// createInteriorVertices
//
// Go through the triangles and create the quad cells
//
// ##################################################################
void createInteriorVertices(GRID * g)
{

   printf("#meshgen: Creating interior vertices ...\n");

   int      i,j,jj,n,nn,is1,is2;
   int      index,halfEdgePerSide,halfVertm1,itemp;
   int      k,kk,k1,k2,k3,k4,k5,k6,ktemp,ktemp3,ktemp4,ktemp5,koffset;
   int      iA,iB,iC,iD,iE,iF,iO;
   int      index1, index2;
   int      edgeID[3],edge[4],dirAB,dirBC,dirCA;
   int     *vertID, iSurface;
   int      nodeMidDO,nodeMidEO,nodeMidFO;
   double   aa,bb,cc;
   double   A[3],B[3],C[3],O[3];
   double   deltaDist;
   double  *boundaryPts, *intPts, *posVert, *normalVert;
   double   dummy1[3],dummy2[3];

   for (i = 0; i < 3; i++)
      dummy1[i] = dummy2[i] = 0.0;

   g->midEdgeID = (int *) malloc(sizeof(int)*3*g->numTriangle);

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
   k6 = 2*g->numTriEdge; // for additional subdivision at triangular level

   aa = 0.5;
   bb = 0.5;
   cc = 0.5;

   n = 3*pow4*g->numTriangle;
   g->numQuadConn = n;

   // ===============================================================
   // Allocations
   // ===============================================================

   g->quadConn     = (int **) malloc(sizeof(int *)*n);
   g->quad2triList = (int *) malloc(sizeof(int)*n);

   for (i = 0; i < n; i++)
      g->quadConn[i]     = (int *) malloc(sizeof(int)*(N_EDGE));
   
   // initialize quadConn
   for (i = 0; i < n; i++)
      for (j = 0; j < N_EDGE; j++)      
         g->quadConn[i][j] = 0;

   boundaryPts = (double *) malloc(sizeof(double)*9*g->nVertPerSide); // 3 sides 
   intPts      = (double *) malloc(sizeof(double)*3*g->nIntPts);

   // allocate posVert   
   deltaDist   = 1./halfEdgePerSide;
   posVert     = (double *) malloc(sizeof(double)*3*halfVertm1);
   normalVert  = (double *) malloc(sizeof(double)*3*halfVertm1);

   vertID      = (int *) malloc(sizeof(int)*(halfVertp2));

   // allocate triQuadEdge
   nn = 2*g->numTriEdge + 3*g->numTriangle; // total number of triangular
                                            // level edges
   g->numTriQuadEdge    = nn;
   g->triQuadEdge       = (int **) malloc(sizeof(int *)*nn);

   for (i = 0; i < nn; i++)
      g->triQuadEdge[i] = (int *) malloc(sizeof(int)*4);

   // initial values
   for (i = 0; i < nn; i++)
      for (j = 0; j < 4; j++)
         g->triQuadEdge[i][j] = -1;

   // ===============================================================
   //
   // Loop over all the triangles and build the cell connectivity
   // and edge information
   //
   // ===============================================================
   for (i = 0; i < g->numTriangle; i++)
   {
      iA = g->triConn[3*i];
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

         // quadEdge: (refer hardcopy documentation)
         // The 3rd and 4th column refers to the index in Qconn array which
         // contains the particular edge
         // The 5th and 6th column refers to the edge number within the said
         // loop - therefore runs from 1-4 (part of a quad cell)
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

               ktemp = is2;
               itemp = 0; // start
               for (k = 0; k < halfEdgePerSide; k++)
               {
                  g->quadEdge[ktemp][3] = k2+itemp;
                  g->quadEdge[ktemp][5] = 0;
                  ktemp++;
                  itemp++; // skip

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]  ][0] = iA;
                  g->triQuadEdge[2*edgeID[j]  ][1] = iD;
                  g->triQuadEdge[2*edgeID[j]  ][3] = 1;

               } // k loop

               itemp = pow4;
               for (k = halfEdgePerSide; k < g->nEdgePerSide; k++)
               {

                  g->quadEdge[ktemp][2] = k2+itemp;
                  g->quadEdge[ktemp][4] = 0;
                  ktemp++;
                  itemp++; // skip

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]+1][0] = iB;
                  g->triQuadEdge[2*edgeID[j]+1][1] = iD;
                  g->triQuadEdge[2*edgeID[j]+1][2] = 1;
               } // k loop

            }
            else if (edge[1] == iC)
            {
               dirCA = -1;
               iF    = is1 + halfVertm1;

               ktemp = is2;
               itemp = 0;
               for (k = 0; k < halfEdgePerSide; k++)
               {
                  g->quadEdge[ktemp][2] = k2+itemp;
                  g->quadEdge[ktemp][4] = 3;
                  ktemp++;
                  itemp += pow2;

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]  ][0] = iA;
                  g->triQuadEdge[2*edgeID[j]  ][1] = iF;
                  g->triQuadEdge[2*edgeID[j]  ][2] = 1;
                  
               } // k loop

               itemp = 2*pow4 + (pow2-1)*pow2;
               for (k = halfEdgePerSide; k < g->nEdgePerSide; k++)
               {

                  g->quadEdge[ktemp][3] = k2+itemp;
                  g->quadEdge[ktemp][5] = 2;
                  ktemp++;
                  itemp++; // skip

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]+1][0] = iC;
                  g->triQuadEdge[2*edgeID[j]+1][1] = iF;
                  g->triQuadEdge[2*edgeID[j]+1][3] = 1;
               } // k loop
               
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
               
               ktemp = is2;
               itemp = pow4 + pow2 - 1;
               for (k = 0; k < halfEdgePerSide; k++)
               {
                  g->quadEdge[ktemp][2] = k2+itemp;
                  g->quadEdge[ktemp][4] = 0;
                  ktemp++;
                  itemp--;

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]  ][0] = iB;
                  g->triQuadEdge[2*edgeID[j]  ][1] = iD;
                  g->triQuadEdge[2*edgeID[j]  ][2] = 1;
               } // k loop

               itemp = pow2 - 1;
               for (k = halfEdgePerSide; k < g->nEdgePerSide; k++)
               {

                  g->quadEdge[ktemp][3] = k2+itemp;
                  g->quadEdge[ktemp][5] = 0;
                  ktemp++;
                  itemp--; // skip

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]+1][0] = iA;
                  g->triQuadEdge[2*edgeID[j]+1][1] = iD;
                  g->triQuadEdge[2*edgeID[j]+1][3] = 1;
               } // k loop

               
            }
            else if (edge[1] == iC)
            {
               dirBC = 1;
               iE    = is1 + halfVertm1;

               ktemp = is2;
               itemp = pow4 + pow2 - 1;
               for (k = 0; k < halfEdgePerSide; k++)
               {
                  g->quadEdge[ktemp][3] = k2+itemp;
                  g->quadEdge[ktemp][5] = 1;
                  ktemp++;
                  itemp += pow2;

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]  ][0] = iB;
                  g->triQuadEdge[2*edgeID[j]  ][1] = iE;
                  g->triQuadEdge[2*edgeID[j]  ][3] = 1;
               } // k loop

               itemp = 2*pow4 + pow2 - 1;
               for (k = halfEdgePerSide; k < g->nEdgePerSide; k++)
               {

                  g->quadEdge[ktemp][2] = k2+itemp;
                  g->quadEdge[ktemp][4] = 1;
                  ktemp++;
                  itemp += pow2; // skip

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]+1][0] = iC;
                  g->triQuadEdge[2*edgeID[j]+1][1] = iE;
                  g->triQuadEdge[2*edgeID[j]+1][2] = 1;
               } // k loop

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
               
               ktemp = is2;
               itemp = 3*pow4 - 1;
               for (k = 0; k < halfEdgePerSide; k++)
               {
                  g->quadEdge[ktemp][3] = k2+itemp;
                  g->quadEdge[ktemp][5] = 2;
                  ktemp++;
                  itemp--;

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]  ][0] = iC;
                  g->triQuadEdge[2*edgeID[j]  ][1] = iF;
                  g->triQuadEdge[2*edgeID[j]  ][3] = 1;

               } // k loop

               itemp = (pow2-1)*pow2;
               for (k = halfEdgePerSide; k < g->nEdgePerSide; k++)
               {

                  g->quadEdge[ktemp][2] = k2+itemp;
                  g->quadEdge[ktemp][4] = 3;
                  ktemp++;
                  itemp -= pow2; // skip

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]+1][0] = iF;
                  g->triQuadEdge[2*edgeID[j]+1][1] = iA;
                  g->triQuadEdge[2*edgeID[j]+1][2] = 1;
               } // k loop

            }
            else if (edge[1] == iB) // edge is C--B
            {
               dirBC = -1;
               iE    = is1 + halfVertm1;

               ktemp = is2;
               itemp = 3*pow4-1;
               for (k = 0; k < halfEdgePerSide; k++)
               {
                  g->quadEdge[ktemp][2] = k2+itemp;
                  g->quadEdge[ktemp][4] = 1;
                  ktemp++;
                  itemp -= pow2;


                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]  ][0] = iC;
                  g->triQuadEdge[2*edgeID[j]  ][1] = iE;
                  g->triQuadEdge[2*edgeID[j]  ][2] = 1;
               
               } // k loop

               itemp = 2*pow4-1;
               for (k = halfEdgePerSide; k < g->nEdgePerSide; k++)
               {

                  g->quadEdge[ktemp][3] = k2+itemp;
                  g->quadEdge[ktemp][5] = 1;
                  ktemp++;
                  itemp -= pow2; // skip

                  // append to triQuadEdge
                  g->triQuadEdge[2*edgeID[j]+1][0] = iB;
                  g->triQuadEdge[2*edgeID[j]+1][1] = iE;
                  g->triQuadEdge[2*edgeID[j]+1][3] = 1;
               } // k loop

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

          g->allNodePos[3*iO + k] = O[k];
          g->nodeNormal[3*iO + k] = ONE_THIRD*(g->nodeNormal[3*iA + k]
                                             + g->nodeNormal[3*iB + k]
                                             + g->nodeNormal[3*iC + k]);
      }

      // flush point O to surface (only if sphere)
      if(strcmp(surfaceType,"sphere")==0 ||
         strcmp(surfaceType,"robin") ==0 )
      {
         moveToBoundary(dummy1,dummy2,&g->allNodePos[3*iO],
            &g->nodeNormal[3*iO],surfaceType);
      }

      ktemp++; // update counter for center position

      // ============================================================
      // add the edges O-D, O-E and O-F to the edges of triangular 
      // level for additional subdivision
      // ============================================================
      g->triQuadEdge[k6  ][0] = iD;
      g->triQuadEdge[k6  ][1] = iO;
      g->triQuadEdge[k6  ][2] = 1;
      g->triQuadEdge[k6  ][3] = 1;

      g->triQuadEdge[k6+1][0] = iE;
      g->triQuadEdge[k6+1][1] = iO;
      g->triQuadEdge[k6+1][2] = 1;
      g->triQuadEdge[k6+1][3] = 1;

      g->triQuadEdge[k6+2][0] = iF;
      g->triQuadEdge[k6+2][1] = iO;
      g->triQuadEdge[k6+2][2] = 1;
      g->triQuadEdge[k6+2][3] = 1;

      k6+=3; // update counter for triQuadEdge
      // ============================================================
      // Loop over each of the three edges of a triangle
      // and add points along the connecting edges
      //
      // Also update the quadEdge
      // ============================================================
      ktemp3 = k3;
      ktemp4 = k4;
      ktemp5 = k5;

      // loop over each edge
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
         // Update the normal values as well  
         for (k = 0; k < halfVertm1; k++)
         {            
            g->allNodePos[3*ktemp4    ] = posVert[3*k    ];
            g->allNodePos[3*ktemp4 + 1] = posVert[3*k + 1];
            g->allNodePos[3*ktemp4 + 2] = posVert[3*k + 2];

            g->nodeNormal[3*ktemp4    ] = normalVert[3*k    ];
            g->nodeNormal[3*ktemp4 + 1] = normalVert[3*k + 1];
            g->nodeNormal[3*ktemp4 + 2] = normalVert[3*k + 2];

            ktemp4++;
         } // k loop

         //vertID = [i1 k1 k1+1 k1+2 ... k1+(nVertPerSide-1) i2]
         vertID[0] = index; vertID[halfVertp2-1] = iO;

         for (k = 0; k < halfVertm1; k++)
            vertID[k+1] = ktemp4 - halfVertm1 + k;


         if      (index == iD) // bordering Quad1 and Quad2
         {
            
            koffset           = g->numTriEdge*g->nEdgePerSide 
                              + i*(3*halfEdgePerSide+3*edgePerQuad);
   
            g->midEdgeID[3*i] = koffset; // edge of DO
            nodeMidDO         = ktemp4 - halfVertm1 + halfhalfVertm1; //
            for (k = 0; k < halfEdgePerSide; k++)
            {
               g->quadEdge[koffset][0] = vertID[k+1];
               g->quadEdge[koffset][1] = vertID[k];               
               g->quadEdge[koffset][2] = k2 + pow4   + k*pow2;
               g->quadEdge[koffset][3] = k2 + pow2-1 + k*pow2; 
               g->quadEdge[koffset][4] = 3;
               g->quadEdge[koffset][5] = 1; 
               koffset++;
               ktemp3++;
            }
         }
         else if (index == iE) // bordering Quad2 and Quad3
         {
            
            koffset             = g->numTriEdge*g->nEdgePerSide 
                                + i*(3*halfEdgePerSide+3*edgePerQuad)
                                +   halfEdgePerSide;


            g->midEdgeID[3*i+1] = koffset; // edge of EO
            nodeMidEO           = ktemp4 - halfVertm1 + halfhalfVertm1;
            for (k = 0; k < halfEdgePerSide; k++)
            {
               g->quadEdge[koffset][0] = vertID[k+1];
               g->quadEdge[koffset][1] = vertID[k];
               g->quadEdge[koffset][2] = k2 + 2*pow4+pow2-1 - k;
               g->quadEdge[koffset][3] = k2 + 2*pow4-1 - k;
               g->quadEdge[koffset][4] = 0;
               g->quadEdge[koffset][5] = 2;
               koffset++;
               ktemp3++;
            }
         }
         else if (index == iF) // bordering Quad1 and Quad3
         {
            koffset             = g->numTriEdge*g->nEdgePerSide 
                                + i*(3*halfEdgePerSide+3*edgePerQuad)
                                + 2*halfEdgePerSide;

            g->midEdgeID[3*i+2] = koffset; // edge of FO
            nodeMidFO           = ktemp4 - halfVertm1 + halfhalfVertm1;
            for (k = 0; k < halfEdgePerSide; k++)
            {
               g->quadEdge[koffset][0] = vertID[k+1];
               g->quadEdge[koffset][1] = vertID[k];
               g->quadEdge[koffset][2] = k2 + pow2*(pow2-1) + k;
               g->quadEdge[koffset][3] = k2 + 2*pow4+pow2*(pow2-1) - k*pow2; 
               g->quadEdge[koffset][4] = 2;
               g->quadEdge[koffset][5] = 3;
               koffset++;
               ktemp3++;
            }
         }
         else
         {
            traces('Something is missing. Stopping. edgeOperations.c');
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
         // obtained. Now obtain the interior points, form the loops
         // and the edges
         // =========================================================
         // Build middle points and find its position
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

               } // k loop
               if(strcmp(surfaceType,"sphere")==0 ||
                  strcmp(surfaceType,"robin")==0)
                  moveToBoundary(dummy1,dummy2,&g->allNodePos[3*ktemp5]
                     ,&g->nodeNormal[3*ktemp5],surfaceType);

               ktemp5++;
            } // kk loop
         } // jj loop

         // =========================================================
         // Build quadConn
         // =========================================================
         ktemp = k2 + j*pow4;
         for (kk = 1; kk < halfVertp2; kk++)
         {
            for (jj = 1; jj < halfVertp2; jj++)
            {
               g->quadConn[ktemp][0] = quadPtID[jj-1][kk-1];
               g->quadConn[ktemp][1] = quadPtID[jj  ][kk-1];
               g->quadConn[ktemp][2] = quadPtID[jj  ][kk  ];
               g->quadConn[ktemp][3] = quadPtID[jj-1][kk  ];

               // the first column of quad2triList is the ID of the
               // quadrilateral and the second column is the ID of the
               // original triangle. The vertices and the area can be
               // computed with this data
               g->quad2triList[ktemp] = i;

               // update
               ktemp++;
            } // jj loop
         } // kk loop
         
         // =========================================================
         // Build quadEdges
         // Number of vert and horz edges: halfVertm1*halfEdgePerSide
         // =========================================================
         // FOR QUAD 1 (horz edge: right->left, vert edge: down->top)
         // =========================================================
         if(j == 0 )
         {
            // Horizontal edges
            for (jj = 0; jj < halfVertm1; jj++)
            {
               for (kk = 0; kk < halfEdgePerSide; kk++)
               {
                  g->quadEdge[ktemp3][1] = quadPtID[jj+1][kk  ];
                  g->quadEdge[ktemp3][0] = quadPtID[jj+1][kk+1];
                  g->quadEdge[ktemp3][3] = k2 + j*pow4 +  jj   + kk*pow2;
                  g->quadEdge[ktemp3][2] = k2 + j*pow4 +  jj+1 + kk*pow2;
                  g->quadEdge[ktemp3][5] = 1;
                  g->quadEdge[ktemp3][4] = 3;
                  ktemp3++;
               }
            }
            // vertical edges
            for (kk = 0; kk < halfVertm1; kk++)
            {
               for (jj = 0; jj < halfEdgePerSide; jj++)
               {
                  g->quadEdge[ktemp3][0] = quadPtID[jj+1][kk+1];
                  g->quadEdge[ktemp3][1] = quadPtID[jj  ][kk+1];
                  g->quadEdge[ktemp3][2] = k2 + j*pow4 + jj + kk*pow2;
                  g->quadEdge[ktemp3][3] = k2 + j*pow4 + jj + (kk+1)*pow2;
                  g->quadEdge[ktemp3][4] = 2;
                  g->quadEdge[ktemp3][5] = 0;
                  ktemp3++;
               }
            }
         }
         // =========================================================
         // FOR QUAD 2 (horz edge: right->left, vert edge: top->down)
         // =========================================================
         else if (j == 1)
         {
            // Horizontal edges
            for (jj = 0; jj < halfVertm1; jj++)
            {
               for (kk = 0; kk < halfEdgePerSide; kk++)
               {
                  g->quadEdge[ktemp3][1] = quadPtID[jj+1][kk  ];
                  g->quadEdge[ktemp3][0] = quadPtID[jj+1][kk+1];
                  g->quadEdge[ktemp3][3] = k2 + j*pow4 +  jj   + kk*pow2;
                  g->quadEdge[ktemp3][2] = k2 + j*pow4 +  jj+1 + kk*pow2;
                  g->quadEdge[ktemp3][5] = 1;
                  g->quadEdge[ktemp3][4] = 3;
                  ktemp3++;
               }
            }
            // vertical edges
            for (kk = 0; kk < halfVertm1; kk++)
            {
               for (jj = 0; jj < halfEdgePerSide; jj++)
               {
                  g->quadEdge[ktemp3][1] = quadPtID[jj+1][kk+1];
                  g->quadEdge[ktemp3][0] = quadPtID[jj  ][kk+1];
                  g->quadEdge[ktemp3][3] = k2 + j*pow4 + jj + kk*pow2;
                  g->quadEdge[ktemp3][2] = k2 + j*pow4 + jj + (kk+1)*pow2;
                  g->quadEdge[ktemp3][5] = 2;
                  g->quadEdge[ktemp3][4] = 0;
                  ktemp3++;
               }
            }
         }
         // =========================================================
         // FOR QUAD 3 (horz edge: left->right, vert edge: top->down)
         // =========================================================
         else if (j == 2)
         {
            // Horizontal edges
            for (jj = 0; jj < halfVertm1; jj++)
            {
               for (kk = 0; kk < halfEdgePerSide; kk++)
               {
                  g->quadEdge[ktemp3][0] = quadPtID[jj+1][kk  ];
                  g->quadEdge[ktemp3][1] = quadPtID[jj+1][kk+1];
                  g->quadEdge[ktemp3][2] = k2 + j*pow4 +  jj   + kk*pow2;
                  g->quadEdge[ktemp3][3] = k2 + j*pow4 +  jj+1 + kk*pow2;
                  g->quadEdge[ktemp3][4] = 1;
                  g->quadEdge[ktemp3][5] = 3;
                  ktemp3++;
               }
            }
            // vertical edges
            for (kk = 0; kk < halfVertm1; kk++)
            {
               for (jj = 0; jj < halfEdgePerSide; jj++)
               {
                  g->quadEdge[ktemp3][1] = quadPtID[jj+1][kk+1];
                  g->quadEdge[ktemp3][0] = quadPtID[jj  ][kk+1];
                  g->quadEdge[ktemp3][3] = k2 + j*pow4 + jj + kk*pow2;
                  g->quadEdge[ktemp3][2] = k2 + j*pow4 + jj + (kk+1)*pow2;
                  g->quadEdge[ktemp3][5] = 2;
                  g->quadEdge[ktemp3][4] = 0;
                  ktemp3++;
               }
            }
            
         }
         else
         {
            traces(Something wrong. Stopping.);
            exit(1);
         }


      } // j loop (loop over each quad)

      k1++;         // for center
      k2 += 3*pow4; // for number of cells in a triangle
      k3  = ktemp3; // edges inside triangle
      k4  = ktemp4; // vertices on int edges
      k5  = ktemp5; // quad int points

	} // i loop (g->numTriangle)

   g->numQuadEdgeT  = g->numQuadEdge; // temp variable used in loopsAndCells
   g->numQuadEdge   = k3;

   for (i = 0; i < g->numQuadEdge; i++)
   {

      if (g->quadEdge[i][2] == -1 || g->quadEdge[i][2] == -5)
      {
         SWAP(g->quadEdge[i][0],g->quadEdge[i][1]);
         SWAP(g->quadEdge[i][2],g->quadEdge[i][3]);
         SWAP(g->quadEdge[i][4],g->quadEdge[i][5]);
      }

   }

}

// ==================================================================
//
// subdomainEdgeBlanking
//
// ==================================================================
void subdomainEdgeBlanking(GRID *g)
{
   int i,j,k;
   int offset;

   for (i = 0; i < g->boundaryCount; i++)
   {
      offset = g->boundaryID[i]*g->nEdgePerSide;
      for (j = 0; j < g->nEdgePerSide; j++)
      {
         g->quadEdge[offset+j][2] = -5;
         g->quadEdge[offset+j][3] = -5;
         g->quadEdge[offset+j][4] = -5;
         g->quadEdge[offset+j][5] = -5;
      }

   }


}

// ##################################################################
// END OF FILE
// ##################################################################
