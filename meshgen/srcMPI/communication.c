// #################################################################
//
// communication.c
//
// Routines that contain MPI based communication routines
// #################################################################
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "globalVariables.h"
#include "communication.h"

// #################################################################
//
// boundaryNodeConnection_MPI
//
// single thread finds out the connecting nodes to each node
// on the sub-domain boundary
//
// performed by only a single thread
//
// Routine extracts the following two arrays:
//     boundaryNodeID  - array containing the boundary nodes in a
//                       given domain
//     boundaryNodeMap - Array containing both the local ID of the
//                       corresponding node in the neighbouring
//                       domain and the domain ID
// #################################################################
void boundaryNodeConnection_MPI(int gridID) 
{


   printf("#meshgen: Finding boundary connections ...\n");

   int    i,j,k,ii,jj,kk,flag,temp;
   int    idx1,idx2,idx3,idx4,idx5,idx6,node1,node2;
   int    displacement[nGrid];
   double dist;   
   int myNodeID,myDomainID,otherNodeID,otherDomainID;
   int idxMyNodeID,idxOtherNodeID;

   // allocate
   allNumNodePos           = (int *) malloc(sizeof(int)*nGrid);
   allNumNodePosCum        = (int *) malloc(sizeof(int)*nGrid);
   allBoundaryNodeCount    = (int *) malloc(sizeof(int)*nGrid);
   allBoundaryNodeCountCum = (int *) malloc(sizeof(int)*nGrid);

   totalBoundaryNodeCount = 0;
   totalNumNodePos = 0;
   g[gridID].boundaryNodeCount = 0;

   g[gridID].ifBoundaryNode = (int *) malloc(sizeof(int)*g[gridID].numNodePos);

   // initialize
   for (j = 0; j < g[gridID].numNodePos; j++)
      g[gridID].ifBoundaryNode[j] = -1;


   // Loop over each grid and identify the nodes that belong 
   // to the edges that lie on the loop
   for (j = 0; j < g[gridID].numQuadEdge; j++)
   {

      // belongs to sub-domain boundary
      if(g[gridID].quadEdge[j][2] == -5 || g[gridID].quadEdge[j][3] == -5 )
      {
         // get the node IDs
         node1 = g[gridID].quadEdge[j][0];
         node2 = g[gridID].quadEdge[j][1];
      
         if (g[gridID].ifBoundaryNode[node1]==-1)
            g[gridID].boundaryNodeCount++;

         if (g[gridID].ifBoundaryNode[node2]==-1)
            g[gridID].boundaryNodeCount++;

         g[gridID].ifBoundaryNode[node1] = 1;
         g[gridID].ifBoundaryNode[node2] = 1;

      } //if loop
   } // j loop

   // allocate for future use
   g[gridID].boundaryNodeID   = (int *) malloc(sizeof(int)*g[gridID].boundaryNodeCount);
   g[gridID].boundaryNodeMap  = (int *) malloc(sizeof(int)*MAX_CONN*2*g[gridID].boundaryNodeCount);
   g[gridID].revMapNumber     = (int *) malloc(sizeof(int)*MAX_CONN*g[gridID].boundaryNodeCount);
   g[gridID].boundaryNodePos  = (double *) malloc(sizeof(double)*3*g[gridID].boundaryNodeCount);

   temp = 2*MAX_CONN*g[gridID].boundaryNodeCount;
   for (j = 0; j < temp; j++)
      g[gridID].boundaryNodeMap[j] = -1;

   temp = MAX_CONN*g[gridID].boundaryNodeCount;
   for (j = 0;  j < temp; j++) 
      g[gridID].revMapNumber[j] = -1;

   k = 0;
   for (j = 0; j < g[gridID].numNodePos; j++)
   {
      if(g[gridID].ifBoundaryNode[j]==1)
      {
         g[gridID].boundaryNodeID[k] = j;
         g[gridID].boundaryNodePos[3*k  ] = g[gridID].allNodePos[3*j  ];
         g[gridID].boundaryNodePos[3*k+1] = g[gridID].allNodePos[3*j+1];
         g[gridID].boundaryNodePos[3*k+2] = g[gridID].allNodePos[3*j+2];
         k++;
      }
   }


   // send the total number of nodes per subdomain to all procs   
   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Allgather(&g[gridID].boundaryNodeCount, 1, MPI_INT, 
      allBoundaryNodeCount, 1, MPI_INT, MPI_COMM_WORLD);

   MPI_Allgather(&g[gridID].numNodePos, 1, MPI_INT, 
      allNumNodePos, 1, MPI_INT, MPI_COMM_WORLD);

   for (j = 0; j < nGrid; j++)
   {
      totalBoundaryNodeCount += allBoundaryNodeCount[j];
      totalNumNodePos        += allNumNodePos[j];
   }

// ==================================================================
// allocation and initializations
// ==================================================================
   allBoundaryNodeID  = (int *) malloc(sizeof(int)*totalBoundaryNodeCount);
   
   allBoundaryNodePos = (double *) malloc(sizeof(double)*3*totalBoundaryNodeCount);

   int *currentMyID = (int *) malloc(sizeof(int)*totalBoundaryNodeCount);
   int *currentOtherID = (int *) malloc(sizeof(int)*totalBoundaryNodeCount);
   for (j =0; j < totalBoundaryNodeCount; j++)
      currentMyID[j] = currentOtherID[j] = 0;

   allBoundaryMapID   = (int **) malloc(sizeof(int *)*totalBoundaryNodeCount);
   for (i = 0; i < totalBoundaryNodeCount; i++)
   {
      temp = 2*MAX_CONN;
      allBoundaryMapID[i] = (int *) malloc(sizeof(int)*temp);      
      for (j = 0; j < temp; j++)
         allBoundaryMapID[i][j] = -1;
   }

   allNodeRevMap = (int **) malloc(sizeof(int *)*totalBoundaryNodeCount);
   for (i = 0; i < totalBoundaryNodeCount; i++)
   {
      allNodeRevMap[i] = (int *) malloc(sizeof(int)*MAX_CONN);
      for (j = 0; j < MAX_CONN; j++)
         allNodeRevMap[i][j] = -1;
   }

// ==================================================================
// Use MPI_Allgatherv to distribute data to all procs
// ==================================================================
   displacement[0] = 0;
   for (j = 1; j < nGrid; j++)
      displacement[j] = displacement[j-1] + allBoundaryNodeCount[j-1];

   // share boundaryNodeID and boundaryNodePos with all procs
   MPI_Allgatherv(g[gridID].boundaryNodeID, g[gridID].boundaryNodeCount, MPI_INT,
      allBoundaryNodeID, allBoundaryNodeCount, displacement, MPI_INT, MPI_COMM_WORLD);

   // Node positions are three times as big as node ID array
   for (j = 0; j < nGrid; j++)
   {
      allBoundaryNodeCount[j] *= 3;
      displacement[j]         *= 3;
   }

   MPI_Allgatherv(g[gridID].boundaryNodePos, 3*g[gridID].boundaryNodeCount, MPI_DOUBLE,
      allBoundaryNodePos, allBoundaryNodeCount, displacement, MPI_DOUBLE, MPI_COMM_WORLD);

   // Node positions are three times as big as node ID array
   for (j = 0; j < nGrid; j++)
   {
      allBoundaryNodeCount[j] /= 3;
      displacement[j]         /= 3;
   }

   // find the cumulative sum of the boundary node count
   allBoundaryNodeCountCum[0] = allNumNodePosCum[0] = 0;
   for (j = 1; j < nGrid; j++)
   {
         allBoundaryNodeCountCum[j] = allBoundaryNodeCountCum[j-1] + allBoundaryNodeCount[j-1];
         allNumNodePosCum[j]        = allNumNodePosCum[j-1] + allNumNodePos[j-1];
   }

   // find the maximum number of nodes of all domains
   int tempint2[nGrid];

   MPI_Gather(&g[gridID].numNodePos, 1, MPI_INT, tempint2, 1, MPI_INT, 0, MPI_COMM_WORLD);

   if (gridID == 0)
   {
      maxNumNodes=0;
      for (i = 0; i < nGrid; i++)
      {
         maxNumNodes = (maxNumNodes > tempint2[i]) ? maxNumNodes : tempint2[i];
      }
   } // if gridID == 0

   MPI_Bcast(&maxNumNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

// ==================================================================
// identify mapping between the different node points
// map the ID's of the node
// ==================================================================   
   MPI_Barrier(MPI_COMM_WORLD);

   allFwdMap = (int **) malloc(sizeof(int *)*nGrid);
   for (i = 0; i < nGrid; i++)
   {
      allFwdMap[i] = (int *) malloc(sizeof(int)*maxNumNodes);
      for (j = 0; j < maxNumNodes; j++)
         allFwdMap[i][j] = -1;
   }
   if (gridID == 0 )
   {
      // create forward MAP, i.e., line number of which row in the
      // 'boundaryNodeMap' a Node belongs to given its ID
      // int **allFwdMap;


      

      for (i = 0; i < nGrid; i++)
      {
         temp = 0;
         for (ii = 0; ii < allBoundaryNodeCount[i]; ii++)
         {
            idxMyNodeID    = allBoundaryNodeCountCum[i] + ii;
            myNodeID       = allBoundaryNodeID[idxMyNodeID];                        

            allFwdMap[i][myNodeID] = temp;
            temp++;
         }
      }


      for (i = 0; i < nGrid; i++)
      {
        printf("Proc 0 performing boundary checks on grid %d\n",i);
        temp = 0;
         // check for absolute distance between two node points
         for (ii = 0; ii < allBoundaryNodeCount[i]; ii++)
         {
            // index of one of the nodes
            idx1 = 3*allBoundaryNodeCountCum[i] + 3*ii;

            flag = 0; // debug flag.. has to change to 1
            for (j = 0; j < nGrid; j++)
            {

               if (i != j) // if not the same grids
               {
                  for (jj = 0; jj < allBoundaryNodeCount[j]; jj++)
                  {
                     
                     idx2 = 3*allBoundaryNodeCountCum[j] + 3*jj;
                     dist = fabs(allBoundaryNodePos[idx1  ]-allBoundaryNodePos[idx2  ])
                          + fabs(allBoundaryNodePos[idx1+1]-allBoundaryNodePos[idx2+1])
                          + fabs(allBoundaryNodePos[idx1+2]-allBoundaryNodePos[idx2+2]);

                     // not the best way to check... :|
                     if(dist < 1.e-14)
                     {
                        idxMyNodeID    = allBoundaryNodeCountCum[i] + ii;
                        idxOtherNodeID = allBoundaryNodeCountCum[j] + jj;

                        myNodeID       = allBoundaryNodeID[idxMyNodeID];                        
                        otherNodeID    = allBoundaryNodeID[idxOtherNodeID];

                        myDomainID     = i;
                        otherDomainID  = j;

                        // idx6 =   allBoundaryNodeCountCum[i] +   ii; 
                        // idx3 = 2*allBoundaryNodeCountCum[i] + 2*ii;
                        idx3 =   2*currentMyID[idxMyNodeID];

                        // idx5 = allBoundaryNodeCountCum[i] + ii;
                        // idx5 = allNumNodePosCum[i] + allBoundaryNodeID[idx5];

                        // idx4 = allBoundaryNodeCountCum[j] + jj;
                        // idx4 = allNumNodePosCum[j] + allBoundaryNodeID[idx4];


                        allBoundaryMapID[idxMyNodeID][idx3  ]  = otherDomainID; //j;  // "other" grid ID
                        allBoundaryMapID[idxMyNodeID][idx3+1]  = otherNodeID; //allBoundaryNodeID[idx4]; // local ID on other grid
                        allNodeRevMap[idxMyNodeID][currentMyID[idxMyNodeID]] =
                           allFwdMap[otherDomainID][otherNodeID];
                        // allNodeRevMap[idxOtherNodeID][currentOtherID[idxOtherNodeID]] = temp;
                        // if(1)
                        // {  
                        //    // trace(idxOtherNodeID);
                        //    // trace(idx3);
                        //    printf("[%d, %d] - > [%d %d || %d %d]\t __ REV(%d):[%d %d]\n",
                        //       myDomainID,myNodeID,
                        //       allBoundaryMapID[idxMyNodeID][0],
                        //       allBoundaryMapID[idxMyNodeID][1],
                        //       allBoundaryMapID[idxMyNodeID][2],
                        //       allBoundaryMapID[idxMyNodeID][3], 
                        //       idxOtherNodeID,
                        //       allNodeRevMap[idxMyNodeID][0],
                        //       allNodeRevMap[idxMyNodeID][1]);                        
                        // }

                        temp++;
                        currentMyID[idxMyNodeID]++;
                        currentOtherID[idxOtherNodeID]++;
                        // flag = 1;
                     } // if dist < 1e-14
                     if(flag==1)temp++;
                  } // jj loop
               } // if loop
            } // j loop

            // if(flag==0)
            // {
            //    printf("Mapping not found for a node --- "
            //       "Grid, node :: (%d, %d)\n",i,ii);
            //    printf("Check this problem. Stopping.\n");
            //    exit(1);
            // }

         } // ii loop
      } // i loop

      
      printf("#meshgen: Proc 0 done finding boundary connections...\n");
   } // if loop (gridID == 0)



// ==================================================================
// Use MPI_scatterv to break data and distribute to the 
// individual procs
// ================================================================== 

   MPI_Barrier(MPI_COMM_WORLD);

   int sendbuf1[2*MAX_CONN*totalBoundaryNodeCount];
   ii = 0;
   for (i = 0; i < totalBoundaryNodeCount; i++)
   {
      for (j = 0; j < 2*MAX_CONN; j++)
      {
         sendbuf1[ii] = allBoundaryMapID[i][j];
         ii++;
      }
   }


   int sendcounts[nGrid];
   for (i = 0; i < nGrid; i++)
      sendcounts[i] = 2*MAX_CONN*allBoundaryNodeCount[i];


   displacement[0] = 0;
   for (j = 1; j < nGrid; j++)
      displacement[j] = displacement[j-1] + sendcounts[j-1];      

// ==================================================================
// Form the boundaryNodeMapID using MPI_Scatterv
// ==================================================================   
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Scatterv(sendbuf1, sendcounts, displacement, MPI_INT,
      g[gridID].boundaryNodeMap, 2*MAX_CONN*g[gridID].boundaryNodeCount,
      MPI_INT, 0, MPI_COMM_WORLD);   


   int sendbuf2[MAX_CONN*totalBoundaryNodeCount];
   ii = 0;
   for (i = 0; i < totalBoundaryNodeCount; i++)
   {
      for (j = 0; j < MAX_CONN; j++)
      {
         sendbuf2[ii] = allNodeRevMap[i][j];
         ii++;
      }
   }

   for (i = 0; i < nGrid; i++)
      sendcounts[i] = MAX_CONN*allBoundaryNodeCount[i];

   displacement[0] = 0;
   for (j = 1; j < nGrid; j++)
      displacement[j] = displacement[j-1] + MAX_CONN*allBoundaryNodeCount[j-1];      

// ==================================================================
// MPI_Scatterv
// ================================================================== 
   MPI_Scatterv(sendbuf2, sendcounts, displacement, MPI_INT,
     g[gridID].revMapNumber, MAX_CONN*g[gridID].boundaryNodeCount, MPI_INT, 0, MPI_COMM_WORLD);   

   MPI_Barrier(MPI_COMM_WORLD);
// ==================================================================
// MPI_Bcast
// broadcast allFwdMap (reforming into single array and then broadcasting it)
// ==================================================================    
   int *tempfwdmap = (int *) malloc(sizeof(int)*nGrid*maxNumNodes);
   k = 0;
   for (i = 0; i < nGrid; i++)
   {
      for (j = 0; j < maxNumNodes; j++)
      {
         tempfwdmap[k] = allFwdMap[i][j];
         k++;
      }
   }

   MPI_Bcast(tempfwdmap, maxNumNodes*nGrid, MPI_INT, 0, MPI_COMM_WORLD);
   // MPI_Bcast(&allFwdMap[0][0], maxNumNodes*nGrid, MPI_INT, 0, MPI_COMM_WORLD);
   k = 0;
   for (i = 0; i < nGrid; i++)
   {
      for (j = 0; j < maxNumNodes; j++)
      {
         allFwdMap[i][j] =tempfwdmap[k];
         k++;
      }
   }
   // if(gridID==0)
   // {
   //    for (j = 0;  j < g[gridID].boundaryNodeCount; j++) 
   //    {            
   //       printf("NODE(%d): boundaryNodeMap[%d %d || %d %d]: revMapNumber[%d %d]\n",j,
   //          g[gridID].boundaryNodeMap[4*j  ],
   //          g[gridID].boundaryNodeMap[4*j+1],
   //          g[gridID].boundaryNodeMap[4*j+2],
   //          g[gridID].boundaryNodeMap[4*j+3],
   //          g[gridID].revMapNumber[2*j  ],
   //          g[gridID].revMapNumber[2*j+1]);
   //    }
   // }

   // ===============================================================
   // SUBDOMAINCONN.DAT
   //
   // Contains the nodes which are connected to each other. Used
   // to prevent artificial discontinuties across subdomains in 
   // tecplot mapping
   // =============================================================== 
   int stridem1,stride;
   int tempnode[MAX_CONN],tempdomain[MAX_CONN];
   int nstrands, otherDomainTemp;
   nstrands = g[gridID].numStrandLayer;

   if (iStrand)
      temp = totalNumNodePos*(nstrands);
   else
      temp = totalNumNodePos;
   
   subdomainconn    = (int *) malloc(sizeof(int)*temp);
   subdomainOtherID = (int *) malloc(sizeof(int)*temp);

   // printf("grid %d, nstrands %d, totalNumNodePos %d, temp %d\n",
   //    gridID,nstrands,totalNumNodePos,temp);

   
   for (i = 0; i < temp; i++)
   {
      subdomainconn[i] = subdomainOtherID[i] =  -1;
   }


   // if(gridID==0)
   // {
   //    printf("gridID %d, node 5, pos: %lf %lf %lf\n",gridID,
   //    g[gridID].allNodePos[3*5],g[gridID].allNodePos[3*5+1],g[gridID].allNodePos[3*5+2] );
   // }
   // if(gridID==1)
   // {
   //    printf("gridID %d, node 32, pos: %lf %lf %lf\n",gridID,
   //    g[gridID].allNodePos[3*32],g[gridID].allNodePos[3*32+1],g[gridID].allNodePos[3*32+2] );
   // }




   if (gridID==0)
   {

      // 
      // ad-hoc fix
      //
      int *tempstride = (int *) malloc(sizeof(int)*nGrid);
      if(iStrand)
      {
         for (i = 0; i < nGrid; i++)
            tempstride[i] = allNumNodePosCum[i]*nstrands;
      }
      else
      {
         for (i = 0; i < nGrid; i++)
            tempstride[i] = allNumNodePosCum[i];  
      }
      //
      // end ad-hoc fix
      //

      ii = 0;
      for (i = 0; i < nGrid; i++)
      {
         // trace(allBoundaryNodeCount[i]);
         // trace(allNumNodePos[i]);
         for (j = 0; j < allBoundaryNodeCount[i]; j++)
         {
            // myNodeID = allNumNodePosCum[i] + allBoundaryNodeID[ii];
            myNodeID = tempstride[i] + allBoundaryNodeID[ii];

            for (k = 0; k < MAX_CONN; k++)
            {
               idx1            = allBoundaryNodeCountCum[i] + j;
               otherDomainTemp = allBoundaryMapID[idx1][2*k];
               otherNodeID     = allBoundaryMapID[idx1][2*k+1];
               if (otherDomainTemp >= 0)
               {
                  tempnode[k] = tempstride[otherDomainTemp] + otherNodeID;
                  // tempnode[k] = allNumNodePosCum[otherDomainTemp] + otherNodeID;
                  otherDomainID = otherDomainTemp;
                  tempdomain[k] = otherDomainID;
               }
               else
               {
                  tempnode[k] = -1;
                  tempdomain[k] = -1;
               }
               
               // printf("myNodeID: %d (%d), otherDomainID: %d, otherNodeID: (%d) %d\n",
               //    myNodeID,i,otherDomainTemp,otherNodeID, tempnode[k]);
               
            }

            // bring the smallest connecting node to the first index
            for (k = MAX_CONN-1; k > 0; k--)
            {
               if (tempnode[k] < 0) continue; // skip the loop

               if (tempnode[k] < tempnode[k-1])
               {
                  temp          = tempnode[k-1];
                  tempnode[k-1] = tempnode[k];
                  tempnode[k]   = temp;

                  temp          = tempdomain[k-1];
                  tempdomain[k-1] = tempdomain[k];
                  tempdomain[k]   = temp;

               }

            }
            if(tempnode[0]==-1)
            {
               traces(something off. stopping);exit(1);
            }
            // printf("myNodeID: %d (%d), otherNodeID: %d (%d)\n",
            //       myNodeID,i,tempnode[0],tempdomain[0]);


            // update subdomainconn array
            if(myNodeID > tempnode[0])
            {
               subdomainconn[myNodeID] = tempnode[0];
               subdomainOtherID[myNodeID] = tempdomain[0];

               if(iStrand)
               {
                  for (jj = 1; jj < nstrands; jj++)
                  {
                     subdomainconn[myNodeID + jj*allNumNodePos[i]] = 
                        tempnode[0] + jj*allNumNodePos[tempdomain[0]];

                     subdomainOtherID[myNodeID + jj*allNumNodePos[i]] = 
                        tempdomain[0];

                  } // jj loop

               } // if loop
            }

            ii++;
         } // j loop
      } // i loop

      //
      // if strands are present
      //
      // if(iStrand)
      // {

      //    for (i = 0; i < totalNumNodePos; i++)
      //    {



      //    }

      //    for (i = 1; i < nstrands-1; i++)
      //    {
      //       // total number of cells in a given layer of the domain
      //       stride   = i*totalNumNodePos;
      //       stridem1 = (i-1)*totalNumNodePos;
      //       for (j = 0; j < totalNumNodePos; j++)
      //       {
      //          if(subdomainconn[j+stridem1] != -1)
      //          {
      //             subdomainconn[j + stride] = subdomainconn[j + stridem1] 
      //                                        + totalNumNodePos;
      //          }
      //       }

      //    }

      // }

      // writearrayINT(subdomainconn,temp);

   }


   // trace(totalNumNodePos);

   // MPI_Abort(MPI_COMM_WORLD,9);

}
// ##################################################################
// END OF FILE
// ##################################################################
