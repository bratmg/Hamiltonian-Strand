// ##################################################################
//
// hybridOperations.c
// 
// Routines for patching the body-fitted structured and outer
// unstructured routines
// ##################################################################

#include "globalVariables.h"
#include "meshtype.h"
#include "hybridOperations.h"

//#define IDX(row,col,ld) (col*l)
#define EPS 0.00001

void meshPatching(HGRID *h)
{
   //
   // variable declarations
   //
   int    i,j,k,ii,jj,kk,stride;
   int    n1,n2,n3,n4,flagn1,flagn2;
   int    idx1,idx2,fourQLevel;
   double x1,y1,z1,x2,y2,z2;

   fourQLevel = 4*qLevel;

// ==================================================================
//
// create the interior points within each segment
//
// ==================================================================
   stride = fourQLevel*(h->nPsi-1)+1;
   for (i = 0; i < h->nEta; i++)
   {
      for (j = 0; j < h->nPsi-1; j++)
      {
         idx1 = stride*i + fourQLevel*j;
         idx2 = stride*i + fourQLevel*(j+1);

         x1   = h->nodePos[3*idx1  ];
         y1   = h->nodePos[3*idx1+1]; 
         z1   = h->nodePos[3*idx1+2];

         x2   = h->nodePos[3*idx2  ];
         y2   = h->nodePos[3*idx2+1]; 
         z2   = h->nodePos[3*idx2+2];

         // number of interior points along an edge
         // depends on the quad-level
         for (k = 1; k < fourQLevel; k++)
         {
            h->nodePos[3*(idx1+k)  ] = x1 + (double)k/fourQLevel*(x2-x1);
            h->nodePos[3*(idx1+k)+1] = y1 + (double)k/fourQLevel*(y2-y1);
            h->nodePos[3*(idx1+k)+2] = z1 + (double)k/fourQLevel*(z2-z1);

         } // k loop

      } // j loop
   } // i loop
   //
   // redefine nPsi
   //
   h->nPsi = (h->nPsi-1)*(fourQLevel) + 1;
   //
   // Link h->nodePos to g[i].allNodePos
   //
   if (nGrid > 1)
   {
      printf("This section had not been extended to multiple nGrid.\n");
      printf("Exiting.\n");
      exit(1);
   }
   // assign nPsi and nEta for later use in GRID *
   for (i = 0; i < nGrid; i++)
   {
      g[i].nPsi = h->nPsi;
      g[i].nEta = h->nEta;
   }
   //
   // append nodes to allNodePos (GRID from HGRID)
   //
   for (i = 0; i < nGrid; i++)
   {
      for (j = 0; j < g[i].numStructNode; j++)
      {
         g[i].allNodePos[3*(g[i].numHamNode+j)  ] = h->nodePos[3*j  ];
         g[i].allNodePos[3*(g[i].numHamNode+j)+1] = h->nodePos[3*j+1];
         g[i].allNodePos[3*(g[i].numHamNode+j)+2] = h->nodePos[3*j+2];

      } // j loop
   } // i loop



// ==================================================================
//
// Find the repeating conn values
//
// ==================================================================
   //
   // isolated edges and nodes that are on boundaries
   //
   for (i = 0; i < nGrid; i++)
   {
      // allocate and initialize
      g[i].numBoundaryEdges = 0;
      g[i].numBoundaryNodes = 0;
      g[i].boundaryEdgeDomain = (int *) malloc(sizeof(int)*g[i].numQuadEdgeHam);
      g[i].boundaryNodeDomain = (int *) malloc(sizeof(int)*8*g[i].numTriNode);
      for (j = 0; j < g[i].numQuadEdgeHam; j++)  g[i].boundaryEdgeDomain[j] = -1;
      for (j = 0; j < 8*g[i].numTriNode; j++)    g[i].boundaryNodeDomain[j] = -1;

      k = 0;
      for (j = 0; j < g[i].numQuadEdgeHam; j++)
      { 
         if(g[i].quadEdge[j][4] == -1 || g[i].quadEdge[j][5] == -1)
         {
            g[i].boundaryEdgeDomain[k] = j;

            // node ID of the edge
            n1 = g[i].quadEdge[j][0];
            n2 = g[i].quadEdge[j][1];

            // check if node already exists
            flagn1 = flagn2 = 0;
            for (ii = 0; ii < g[i].numBoundaryNodes; ii++)
            {
               // node1 already exits
               if (g[i].boundaryNodeDomain[ii] == n1)
                  flagn1 = 1;

               // node2 already exits
               if (g[i].boundaryNodeDomain[ii] == n2)
                  flagn2 = 1;

               // if noth nodes already exist
               if (flagn1 == 1 && flagn2 == 1)
                  break;

            } // ii loop

            // node n1 does not exist
            if (flagn1 == 0)
            {
               g[i].boundaryNodeDomain[g[i].numBoundaryNodes] = n1;
               g[i].numBoundaryNodes++;
            }

            // node n2 does not exist
            if (flagn2 == 0)
            {
               g[i].boundaryNodeDomain[g[i].numBoundaryNodes] = n2;
               g[i].numBoundaryNodes++;
            }

            k++;            
         }
      } // j loop

      g[i].numBoundaryEdges = k;

   } // i loop

   //
   // Compare the nodes with those in boundaryNodeDomain
   // and find the corresponding grid and conn values
   //
   int    flag;
   double dist,distMin,xHam,yHam,zHam;

   // nodesub = [ DOMAIN ID | NODE ID ... ]
   h->nodesub = (int *) malloc(sizeof(int)*h->nPsi);
   for (i = 0; i < h->nPsi; i++) h->nodesub[i] = -1;

   for (j = 0; j < h->nPsi; j++)
   {
      idx1 = h->nPsi*(h->nEta-1) + j;
      x1 = h->nodePos[3*idx1  ];
      y1 = h->nodePos[3*idx1+1];
      z1 = h->nodePos[3*idx1+2];

      distMin = 1.e10;
      flag    = 1;

      //
      // Loop through the grids trying to match nodes
      //
      for (i = 0; i < nGrid; i++)
      {  
         for (ii = 0; ii < g[i].numBoundaryNodes; ii++)
         {
            xHam = g[i].allNodePos[3*g[i].boundaryNodeDomain[ii]  ];
            yHam = g[i].allNodePos[3*g[i].boundaryNodeDomain[ii]+1];
            zHam = g[i].allNodePos[3*g[i].boundaryNodeDomain[ii]+2];

            dist = fabs(xHam-x1) + fabs(yHam-y1) + fabs(zHam-z1);

            if (dist < distMin) distMin = dist;

            if (dist < EPS)
            {  
               // h->nodesub[2*j  ] = i;
               // h->nodesub[2*j+1] = g[i].boundaryNodeDomain[ii];
               h->nodesub[j] = g[i].boundaryNodeDomain[ii];
               flag              = 0;  
            }
         } // ii loop
         
         // error check
         if (flag == 1)
         {            
            printf("Node not linked between interior and exterior mesh\n");
            printf("Stopping. Check error.\n");
            exit(1);
         }

      } // i loop
   } // j loop
   //
   // Update connectivity
   //
   if (nGrid > 1)
   {
      printf("This section has not been updated for multiple nGrid files\n");
      printf("Exiting.\n");
      exit(1);
   }
   int count,offset;
   int leftidx,rightidx,flip;
   int n1_struct,n2_struct,n1_ham,n2_ham,edge_struct,edge_ham;

   int   tempint;
   int * idx1_all = (int *) malloc(sizeof(int)*h->nPsi);
   int * idx2_all = (int *) malloc(sizeof(int)*h->nPsi);
   int * idx3_all = (int *) malloc(sizeof(int)*h->nPsi);
   int * temp_idx = (int *) malloc(sizeof(int)*h->nPsi);


   h->edgesub = (int *) malloc(sizeof(int)*(h->nPsi-1));

   // intialize
   for (j = 0; j < h->nPsi-1; j++)
      h->edgesub[j] = -1;

   for (i = 0; i < nGrid; i++)
   {
      //
      // CREATE h->quadConn ARRAY
      //      

      // define stride lengths
      h->nodeStride = g[i].numHamNode;// g->numNodePos;
      h->connStride = g[i].numQuadConnHam; //g[i].numQuadConn;
      h->edgeStride = g[i].numQuadEdgeHam; //g->numQuadEdge;

      stride = (h->nPsi-1)*(h->nEta-1);

      h->quadConn     = (int **) malloc(sizeof(int *)*stride);

      for (j = 0; j < stride; j++)
         h->quadConn[j] = (int *) malloc(sizeof(int)*(N_EDGE));
      
      // initialize quadConn
      for (j = 0; j < stride; j++)
         for (k = 0; k < N_EDGE; k++)      
            h->quadConn[j][k] = 0;

      //
      // create the h->quadConn arrays
      //
      kk = 0;
      for (j = 0; j < h->nEta-1; j++)
      {
         for (k = 0; k < h->nPsi-1; k++)
         {
            n1 = (j  )*h->nPsi +  k    + h->nodeStride;
            n2 = (j  )*h->nPsi + (k+1) + h->nodeStride;
            n3 = (j+1)*h->nPsi + (k+1) + h->nodeStride;
            n4 = (j+1)*h->nPsi +  k    + h->nodeStride;

            // stitch between the first and last column
            // (periodic conditions on vertical lines n2->n1 and n3->n4)
            if (k == (h->nPsi-2) && iHybPeriodic)
            {
               n2 = (j  )*h->nPsi + h->nodeStride;
               n3 = (j+1)*h->nPsi + h->nodeStride;
            }

            // stitch the domain between the inner mesh 
            // and the outer unstructured mesh 
            if (j == (h->nEta-2))
            {
               n3 = h->nodesub[k+1];
               n4 = h->nodesub[k  ];
            }

            // build quadconn
            h->quadConn[kk][0] = n1;
            h->quadConn[kk][1] = n2;
            h->quadConn[kk][2] = n3;
            h->quadConn[kk][3] = n4;

            kk++;
         } // k loop
      } // j loop

      h->numQuadConn = kk;

      //
      // CREATE quadEdge array for the structured mesh
      //
      // Stitch the node values as done previously for conn values
      //      
      h->numQuadEdge = h->nEta*(h->nPsi-1) + (h->nEta-1)*h->nPsi;

      // allocate and initialize
      h->quadEdge = (int **) malloc(sizeof(int *)*h->numQuadEdge);
      for (j = 0; j < h->numQuadEdge; j++)
      {
         h->quadEdge[j] = (int *) malloc(sizeof(int)*N_QEDGE);
         for (k = 0; k < N_QEDGE; k++)
            h->quadEdge[j][k] = -1;
      }
      // 
      // create horizontal edge arrays
      //
      kk = 0;
      for (j = 0; j < h->nEta; j++)
      {
         for (k = 0; k < h->nPsi-1; k++)
         {
            // left and right node
            h->quadEdge[kk][0] = j*h->nPsi + k   + h->nodeStride;
            h->quadEdge[kk][1] = j*h->nPsi + k+1 + h->nodeStride;

            // periodic boundary conditions for node == h->nPsi
            if (k == (h->nPsi-2) && iHybPeriodic)
            {
               h->quadEdge[kk][1] = j*h->nPsi + h->nodeStride;
            }

            // blending (top edge, match with boundary node)
            if (j == (h->nEta-1))
            {
               h->quadEdge[kk][0] = h->nodesub[k  ];
               h->quadEdge[kk][1] = h->nodesub[k+1];
            }

            // left and right cell
            h->quadEdge[kk][2] = kk               + h->connStride;
            h->quadEdge[kk][3] = kk - (h->nPsi-1) + h->connStride;

            // left and right cell index
            h->quadEdge[kk][4] = 0;
            h->quadEdge[kk][5] = 2;

            // bottom most edge (no right cell)
            if (j == 0)
            {
               h->quadEdge[kk][3] = -1;
               h->quadEdge[kk][5] = -1;
            }

            // top most edge (no left cell)
            if (j == (h->nEta-1))
            {
               h->quadEdge[kk][2] = -1;
               h->quadEdge[kk][4] = -1;
            }

            kk++;

         } // k loop
      } // j loop
      //
      // create vertical edge arrays
      //
      count = 0;
      for (j = 0; j < h->nEta-1; j++)
      {
         count--;
         for (k = 0; k < h->nPsi; k++)
         {

            // left and right node of the edge
            h->quadEdge[kk][1] = (j  )*h->nPsi + k + h->nodeStride;
            h->quadEdge[kk][0] = (j+1)*h->nPsi + k + h->nodeStride;

            // blending for the last node == h->nPsi
            if (k == (h->nPsi-1) && iHybPeriodic)
            {
               h->quadEdge[kk][1] = (j  )*h->nPsi + h->nodeStride;
               h->quadEdge[kk][0] = (j+1)*h->nPsi + h->nodeStride;
            }

            // blending for the top most edge
            if ( j == (h->nEta-2))
            {
               h->quadEdge[kk][0] = h->nodesub[k];
            }

            // left and right cell
            h->quadEdge[kk][3] = count     + h->connStride;
            h->quadEdge[kk][2] = count + 1 + h->connStride;

            // left and right cell index
            h->quadEdge[kk][5] = 1;
            h->quadEdge[kk][4] = 3;

            if(iHybPeriodic == 0)
            {
               // left most edge (no right cell)
               // (recall edge from top to bottom)
               if (k == 0)
               {
                  h->quadEdge[kk][3] = -1;
                  h->quadEdge[kk][5] = -1;
               }

               // right most edge (no left cell)
               // (recall edge from top to bottom)
               if (k == h->nPsi-1)
               {
                  h->quadEdge[kk][2] = -1;
                  h->quadEdge[kk][4] = -1;
               }

            } // not a periodic solution


            count++;
            kk++;

         } // k loop
      } // j loop
      //
      // Map the horizontal edges on the "h->nEta" line
      // with that of the ham mesh boundaries
      //
      for (j = 0; j < h->nPsi-1; j++)
         h->edgesub[j] = -1;

      for (j = 0; j < g[i].numBoundaryEdges; j++)
      {
         n1 = g[i].quadEdge[g[i].boundaryEdgeDomain[j]][0];         
         n2 = g[i].quadEdge[g[i].boundaryEdgeDomain[j]][1];   

         // initialize
         for (k = 0; k < h->nPsi; k++)
            idx1_all[k] = idx2_all[k] = idx3_all[k] = -1;

         findValue(h->nodesub,h->nPsi,n1,idx1_all);
         findValue(h->nodesub,h->nPsi,n2,idx2_all);

         // ad-hoc solution to identify boundary
         if(idx1_all[0] == 2 || idx2_all[0]==2)
         {  
            if (idx1_all[0] > 1 || idx2_all[0] > 2)
            {
               printf("Check this section of the algorithm.\n");
               printf("Stopping\n");
               exit(1);
            }
            idx1 = idx1_all[1];

            for (k = 1; k < h->nPsi; k++)
               temp_idx[k-1] = abs(idx2_all[k]-idx1);

            findValue(temp_idx,h->nPsi-1,1,idx3_all);
            idx2 = idx2_all[idx3_all[1]+1];

         }
         else
         {
            idx1 = idx1_all[1];
            idx2 = idx2_all[1];
         }

         // swap to ensure idx1 < idx2
         if (idx1 > idx2)
         {
            SWAP(idx1,idx2);
            tempint = idx2_all[0];
         }
         else
         {
            tempint = idx1_all[0];
         }

         // h->edgesub runs from 1->(h->nPsi-1) along "h->nEta" line
         
         if (tempint != 0)
         {
            h->edgesub[idx1] = g[i].boundaryEdgeDomain[j];      
         }

      } // j loop

      //
      // error check for h->edgesub
      //
      findValue(h->edgesub,h->nPsi-1,-1,idx3_all);
      if(idx3_all[0] != 0)
      {
         printf("Edge mapping is not accurate. Check it. Stopping.\n");
         exit(1);
      }
      //
      // blend the edges
      //  - Vertical edges: periodic boundary condition (only if iHybPeriodic == 1)
      //
      if(iHybPeriodic)
      {
         offset = h->nEta*(h->nPsi-1);
         for (j = 0; j < h->nEta-1; j++)
         {
            rightidx = j*h->nPsi + h->nPsi-1 + offset;
            leftidx  = j*h->nPsi +             offset;

            // match left cell and right cell
            h->quadEdge[rightidx][2] = h->quadEdge[ leftidx][2];
            h->quadEdge[rightidx][4] = h->quadEdge[ leftidx][4];
            h->quadEdge[ leftidx][3] = h->quadEdge[rightidx][3];
            h->quadEdge[ leftidx][5] = h->quadEdge[rightidx][5];

         } // j loop
      } // if iHybPeriodic
      //
      // blend the edges
      //  - Horizontal edges: the topmost edges (i.e., at h->nEta-1)
      //    are the same as the interior nodes of the ham mesh
      //
      offset = (h->nEta-1)*(h->nPsi-1);
      for (j = 0; j < h->nPsi-1; j++)
      {
         edge_struct = j + offset; // structured mesh edge index
         edge_ham    = h->edgesub[j]; // ham mesh edge index

         // node ID in structured mesh
         n1_struct   = h->quadEdge[edge_struct][0];
         n2_struct   = h->quadEdge[edge_struct][1];      

         // node ID in ham mesh
         n1_ham      = g[i].quadEdge[edge_ham][0];
         n2_ham      = g[i].quadEdge[edge_ham][1];

         dist = fabs(g[i].allNodePos[3*n1_ham  ]-g[i].allNodePos[3*n1_struct  ]) + 
                fabs(g[i].allNodePos[3*n1_ham+1]-g[i].allNodePos[3*n1_struct+1]) + 
                fabs(g[i].allNodePos[3*n1_ham+2]-g[i].allNodePos[3*n1_struct+2]);

         flip = 0;
         if (dist > EPS)
         {
            flip = 1;
            dist = fabs(g[i].allNodePos[3*n1_ham  ]-g[i].allNodePos[3*n2_struct  ]) + 
                   fabs(g[i].allNodePos[3*n1_ham+1]-g[i].allNodePos[3*n2_struct+1]) + 
                   fabs(g[i].allNodePos[3*n1_ham+2]-g[i].allNodePos[3*n2_struct+2]);
            if(dist > EPS)
            {
               trace(j);
               tracef(EPS);
               printf("Nodes do not match. Check. Stopping.\n");
               exit(1);
            }
         }

         //
         // swap the left and right cells for all edges in that column
         //
         if (flip == 1)
         {
            for (k = 0; k < h->nEta; k++)
            {
               idx1 = k*(h->nPsi-1) + j;

               SWAP(h->quadEdge[idx1][0],h->quadEdge[idx1][1]);
               SWAP(h->quadEdge[idx1][2],h->quadEdge[idx1][3]);
               SWAP(h->quadEdge[idx1][4],h->quadEdge[idx1][5]);

            } // k loop
         }

         //
         // blend the edges, match left and right cell
         //
         if (g[i].quadEdge[edge_ham][2]==-1)
         {
            g[i].quadEdge[edge_ham][2]  = h->quadEdge[edge_struct][2];
            g[i].quadEdge[edge_ham][4]  = h->quadEdge[edge_struct][4];
            h->quadEdge[edge_struct][3] = g[i].quadEdge[edge_ham][3];
            h->quadEdge[edge_struct][5] = g[i].quadEdge[edge_ham][5];
         }
         else if (g[i].quadEdge[edge_ham][3]==-1)
         {
            g[i].quadEdge[edge_ham][3]  = h->quadEdge[edge_struct][3];
            g[i].quadEdge[edge_ham][5]  = h->quadEdge[edge_struct][5];
            h->quadEdge[edge_struct][2] = g[i].quadEdge[edge_ham][2];
            h->quadEdge[edge_struct][4] = g[i].quadEdge[edge_ham][4];  
         }

      } // j loop
      // 
      // Append h->quadEdge to g[i].quadEdge array
      //      
      if (nGrid > 1)
      {
         printf("This section had not been extended to multiple nGrid.\n");
         printf("Exiting.\n");
         exit(1);
      }

      // append quadedge
      for (j = 0; j < g[i].numQuadEdgeStruct; j++)
         for (k = 0; k < N_QEDGE; k++)
            g[i].quadEdge[g[i].numQuadEdgeHam+j][k] = h->quadEdge[j][k];

      // append quadconn   
      for (j = 0; j < g[i].numQuadConnStruct; j++)
         for (k = 0; k < N_EDGE; k++)
            g[i].quadConn[g[i].numQuadConnHam+j][k] = h->quadConn[j][k];

      for (j = 0; j < g[i].numQuadEdge; j++)
      {

         if (g[i].quadEdge[j][2] == -1 || g[i].quadEdge[j][2] == -5)
         {
            SWAP(g[i].quadEdge[j][0],g[i].quadEdge[j][1]);
            SWAP(g[i].quadEdge[j][2],g[i].quadEdge[j][3]);
            SWAP(g[i].quadEdge[j][4],g[i].quadEdge[j][5]);
         }

      }

// ==================================================================
//
// Create the extended ham loops on the boundary
//
// ==================================================================
      //
      // Isolate boundary tri nodes
      //
      double posx, posy;
      if(strcmp(surfaceType,"naca") != 0 && 0)
      {
         
         // exit(1);
      }
      else
      {
         g[i].boundaryTriNode   = (int *) malloc(sizeof(int)*g[i].numTriNode);
         g[i].isBoundaryTriNode = (int *) malloc(sizeof(int)*g[i].numTriNode);
         
         for (j = 0; j < g[i].numTriNode; j++)
            g[i].boundaryTriNode[j] = g[i].isBoundaryTriNode[j] = 0;

         k = 0;
         for (j = 0; j < g[i].numBoundaryNodes; j++)
         {
            if (g[i].boundaryNodeDomain[j] < g[i].numTriNode)
            {
               // obtain position
               posx = g[i].allNodePos[3*g[i].boundaryNodeDomain[j]  ];
               posy = g[i].allNodePos[3*g[i].boundaryNodeDomain[j]+1];

               //if ( (fabs(posy) < .62) && (fabs(posx) < 1.1) )
               if ( (fabs(posy) < 2.) && ((posx) < 15.) )
               {
                  g[i].boundaryTriNode[k] = j;
                  g[i].isBoundaryTriNode[g[i].boundaryNodeDomain[j]] = 1;
                  k++;
               }
            }
         } // j loop
         g[i].numBoundaryTriNode = k;
      } // strcmp if loop


      printf("\n\n CHECK BOUNDARY NODE LIMITS! - %d B-NODES\n\n\n",
         g[i].numBoundaryTriNode);


// ==================================================================
//
// Reposition the nodes in the structured mesh to ensure that the
// points on the surface are accurately flushed. Care should be taken
// to reposition all the interal nodes along a "h->neta" line
//
// ==================================================================

      double dummy[3];
      double xold,yold,xnew,ynew,x0,y0;
      double sx,sy; // stretching ratios

      for (j = 0; j < h->nPsi; j++)
      {
         if(strcmp(surfaceType,"naca")==0)
         {

            // ensure it is not on the surface to begin with
            if(j % (4*qLevel) != 0)
            {
               n1 = 4*qLevel*floor((double)j/(double)(4.*qLevel));
               n2 = 4*qLevel*ceil((double)j/(double)(4.*qLevel));

               n1 += h->nodeStride;
               n2 += h->nodeStride;

               // node positions at the boundary/surface end
               tempint = j + h->nodeStride;
               xold    = g[i].allNodePos[3*tempint  ];
               yold    = g[i].allNodePos[3*tempint+1];

               // move the node to the boundary
               moveToBoundary(&g[i].allNodePos[3*n1],
                              &g[i].allNodePos[3*n2], 
                              &g[i].allNodePos[3*(j+h->nodeStride)],
                              dummy,
                              surfaceType);

               // node positions at the boundary/surface end
               tempint = j + h->nodeStride;
               xnew    = g[i].allNodePos[3*tempint  ];
               ynew    = g[i].allNodePos[3*tempint+1];
               
               // node position at the boundary between Ham and struct
               tempint = h->nPsi*(h->nEta-1) + j + h->nodeStride;
               x0      = g[i].allNodePos[3*tempint  ];
               y0      = g[i].allNodePos[3*tempint+1];

               //
               // "Stretch" the other points along the "time" line
               //
               sx = (xnew - x0)/(xold-x0);
               sy = (ynew - y0)/(yold-y0);

               if(sx != sx || sy != sy)
               {
                  printf("\nError in stretching definition.\n");
                  printf("Check division by zero. Stopping.\n");
                  exit(1);
               }

               // reposition the nodes
               for (k = 1; k < h->nEta-1; k++)
               {
                  tempint = h->nPsi*k + j + h->nodeStride;
                  
                  g[i].allNodePos[3*tempint  ] = x0 
                     + sx*(g[i].allNodePos[3*tempint  ] - x0);
                  
                  g[i].allNodePos[3*tempint+1] = y0 
                     + sy*(g[i].allNodePos[3*tempint+1] - y0);

               } // k loop

            } // if loop

         } // if (strcmp naca) loop


      } // j loop

   } // i loop

}

// ##################################################################
//
// findValue
//
// finds the array indices of an array matching a particular value
// ##################################################################
void findValue(const int *array, const int n, const int val, int *nodeID)
{
   int i,count;

   count = 0;
   for (i = 0; i < n; i++)
   {  
      // val exists in array
      if (val == array[i])
      {
         nodeID[count+1] = i;
         count++;
      }
   }

   // first entry in nodeID is the total number of matching vals
   nodeID[0] = count;

}





// ##################################################################
// END OF FILE
// ##################################################################
