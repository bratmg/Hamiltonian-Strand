// ##################################################################
//
// strand.c
// 
// Routines for strand grids
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "globalVariables.h"
#include "meshtype.h"
#include "strand.h"

#define N_CELL_FACE 8
// #################################################################
//
// createStrands
//
// Wrapper function to call other routines
// #################################################################
void createStrands(GRID *g)
{
   //
   createStrandTemplate(g);
   //
   updateVariables(g);
   
}

// #################################################################
//
// createStrandTemplate
//
// strands are defined by their spacing along the strand, 
// the starting point, clipping index, and normal direction
// #################################################################
void createStrandTemplate(GRID *g)
{
   int      i,j,temp;
   double   normDist,stepsize;
   double **normal;

   // create strands 
   g->templateDist = (double *) malloc(sizeof(double)*g->numStrandLayer);

   // grow strands from base (based on init mesh length and growth ratio)
   g->templateDist[0] = 0.;

   for (i = 1; i < g->numStrandLayer; i++)   
   {
      g->templateDist[i] = g->templateDist[i-1] 
                         + g->initMeshLen*pow(g->meshGrowth,i-1);
   }

   //(the number of strands is equal to number of node points)
   g->strandGrid = (STRAND *) malloc(sizeof(STRAND)*g->numNodePos);

   // copy allNodePos to surfNodePos
   temp = 3*g->numNodePos;
   for (i = 0; i < temp; i++) g->surfNodePos[i] = g->allNodePos[i];

   // set the normal vectors based on the normal direction of the point
   // computeNormals(g);

   printf("WARNING!!!!!!!!!!!!!!!!! AD HOC SOLUTION....\n");
   printf("COMPUTE NORMALS HAS BEEN BYPASSED\n");
   for (i = 0; i < g->numNodePos; i++)
   {
      g->strandGrid[i].normal[0] = 0.;
      g->strandGrid[i].normal[1] = 0.;
      g->strandGrid[i].normal[2] = 1.;
   }

   for (i = 0; i < g->numNodePos; i++)
      g->strandGrid[i].pos = (double *) malloc(sizeof(double)*3*g->numStrandLayer);


   // create the grid positions for the strands
   for (i = 0; i < g->numNodePos; i++)
   {
      for (j = 0; j < g->numStrandLayer; j++)
      {

         g->strandGrid[i].pos[3*j  ] = g->allNodePos[3*i  ] 
                                     + g->strandGrid[i].normal[0]*g->templateDist[j];
         g->strandGrid[i].pos[3*j+1] = g->allNodePos[3*i+1] 
                                     + g->strandGrid[i].normal[1]*g->templateDist[j];
         g->strandGrid[i].pos[3*j+2] = g->allNodePos[3*i+2] 
                                     + g->strandGrid[i].normal[2]*g->templateDist[j];

         if(g->strandGrid[i].pos[3*j] != g->strandGrid[i].pos[3*j])
         {
            printf("\n Problems at the strands.. Check this. Stopping\n");
            exit(1);
         }

      } // j loop

   } // i loop




   // ===============================================================
   // New method to create the strands using an "advancing front"
   // ===============================================================
   if(strcmp(surfaceType,"robin")==0)
   {
      printf("#meshgen: WARNING!! ADVANCING FRONT ACTIVE FOR ROBIN.\n");

      g->flagalter = (int *) malloc(sizeof(int)*3*g->numNodePos);

      for (i = 0; i < g->numNodePos; i++)
      {
         g->strandGrid[i].pos[0] = g->surfNodePos[3*i  ];
         g->strandGrid[i].pos[1] = g->surfNodePos[3*i+1];
         g->strandGrid[i].pos[2] = g->surfNodePos[3*i+2];

         // adhoc solution to create the volume mesh
         if (g->strandGrid[i].pos[0] > 0.8 && g->strandGrid[i].pos[0] < 1.8)
            g->flagalter[i] = 1;
         else
            g->flagalter[i] = 0;
      }

      for (j = 1; j < g->numStrandLayer; j++)
      {
         stepsize = g->initMeshLen*pow(g->meshGrowth,j);

         // compute the normals
         averageNormals(g);

         for (i = 0; i < g->numNodePos; i++)
         {
            g->strandGrid[i].pos[3*j  ] = g->surfNodePos[3*i  ]
                                        + g->strandGrid[i].normal[0]*stepsize;
            g->strandGrid[i].pos[3*j+1] = g->surfNodePos[3*i+1]
                                        + g->strandGrid[i].normal[1]*stepsize;
            g->strandGrid[i].pos[3*j+2] = g->surfNodePos[3*i+2]
                                        + g->strandGrid[i].normal[2]*stepsize;

            // update surfNodePos
            g->surfNodePos[3*i  ] = g->strandGrid[i].pos[3*j  ];
            g->surfNodePos[3*i+1] = g->strandGrid[i].pos[3*j+1];
            g->surfNodePos[3*i+2] = g->strandGrid[i].pos[3*j+2];

         } // i loop

      } // j loop

   }


}
// ##################################################################
//
// computeNormals
//
// Compute the normals of the strands grids obtain the normals at
// the different quad panels and at the vertex nodes by averaging
// the surrounding panels
// #################################################################
void computeNormals(GRID *g)
{

   int     i,j,temp;
   int     node1,node2,node3,node4;   
   int    *nodeCount;      
   double  tang1[3],tang2[3],norm[3],dist;      
   double *normal;
      
   nodeCount  = (int *)    malloc(sizeof(int)*g->numNodePos);
   normal     = (double *) malloc(sizeof(double)*3*g->numNodePos);   

   // set all iflag and nodeCount for all nodes to 0
   for (i = 0; i < g->numNodePos; i++) nodeCount[i]=0;

   // loop through all the faces (i.e., edges) to obtain the number
   // of edges connecting any given node
   for (i = 0; i < g->numQuadEdge; i++)
   {
      // two node IDs of a given edge
      node1 = g->quadEdge[i][0];
      node2 = g->quadEdge[i][1];

      nodeCount[node1]++; 
      nodeCount[node2]++;      

   }

   // initialize normal array
   temp = 3*g->numNodePos;
   for (i = 0; i < temp; i++) normal[i]=0.0;
  
   for (i = 0; i < g->numQuadConn; i++)
   {
      node1     = g->quadConn[i][0]; node2 = g->quadConn[i][1];
      node3     = g->quadConn[i][2]; node4 = g->quadConn[i][3];

      // compute the normal of the panel
      tang1[0]  = g->surfNodePos[3*node3  ] - g->surfNodePos[3*node1  ];
      tang1[1]  = g->surfNodePos[3*node3+1] - g->surfNodePos[3*node1+1];
      tang1[2]  = g->surfNodePos[3*node3+2] - g->surfNodePos[3*node1+2];

      tang2[0]  = g->surfNodePos[3*node4  ] - g->surfNodePos[3*node2  ];
      tang2[1]  = g->surfNodePos[3*node4+1] - g->surfNodePos[3*node2+1];
      tang2[2]  = g->surfNodePos[3*node4+2] - g->surfNodePos[3*node2+2];

      norm[0]   =  tang1[1]*tang2[2] - tang1[2]*tang2[1]; //  t1y*t2z - t1z*t2y
      norm[1]   = -tang1[0]*tang2[2] + tang1[2]*tang2[0]; // -t1x*t2z + t1z*t2x
      norm[2]   =  tang1[0]*tang2[1] - tang1[1]*tang2[0]; //  t1x*t2y - t1y*t2x

      // assign this normal to the surrounding nodes
      normal[3*node1  ] += norm[0];
      normal[3*node1+1] += norm[1];
      normal[3*node1+2] += norm[2];

      normal[3*node2  ] += norm[0];
      normal[3*node2+1] += norm[1];
      normal[3*node2+2] += norm[2];

      normal[3*node3  ] += norm[0];
      normal[3*node3+1] += norm[1];
      normal[3*node3+2] += norm[2];

      normal[3*node4  ] += norm[0];
      normal[3*node4+1] += norm[1];
      normal[3*node4+2] += norm[2];

   } // i loop

   // average the normals and assign to strand grids
   for (i = 0; i < g->numNodePos; i++)
   {

      if(nodeCount[i] == 0)
      {
         printf("Something amiss here. nodeCount == 0.\n");
         exit(1);
      }

      g->strandGrid[i].normal[0] = normal[3*i  ]/nodeCount[i];
      g->strandGrid[i].normal[1] = normal[3*i+1]/nodeCount[i];
      g->strandGrid[i].normal[2] = normal[3*i+2]/nodeCount[i];

      // normalize 
      dist = g->strandGrid[i].normal[0]*g->strandGrid[i].normal[0] +
             g->strandGrid[i].normal[1]*g->strandGrid[i].normal[1] +
             g->strandGrid[i].normal[2]*g->strandGrid[i].normal[2];

      dist = 1./sqrt(dist);

      g->strandGrid[i].normal[0]*=dist;
      g->strandGrid[i].normal[1]*=dist;
      g->strandGrid[i].normal[2]*=dist;


   }

   // NOTE: allNodePos not changed to surfNodePos for the case
   //       of a sphere. Not the most cleanest way to write code.
   if(strcmp(surfaceType,"sphere")==0)
   {
      for (i = 0; i < g->numQuadConn; i++)
      {
         node1     = g->quadConn[i][0]; node2 = g->quadConn[i][1];
         node3     = g->quadConn[i][2]; node4 = g->quadConn[i][3];

         norm[0]   = g->allNodePos[3*node1  ];
         norm[1]   = g->allNodePos[3*node1+1];
         norm[2]   = g->allNodePos[3*node1+2];

         dist = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];
         dist = 1./sqrt(dist);

         g->strandGrid[node1].normal[0] = dist*norm[0];
         g->strandGrid[node1].normal[1] = dist*norm[1];
         g->strandGrid[node1].normal[2] = dist*norm[2];

         norm[0]   = g->allNodePos[3*node2  ];
         norm[1]   = g->allNodePos[3*node2+1];
         norm[2]   = g->allNodePos[3*node2+2];

         dist = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];
         dist = 1./sqrt(dist);

         g->strandGrid[node2].normal[0] = dist*norm[0];
         g->strandGrid[node2].normal[1] = dist*norm[1];
         g->strandGrid[node2].normal[2] = dist*norm[2];

         norm[0]   = g->allNodePos[3*node3  ];
         norm[1]   = g->allNodePos[3*node3+1];
         norm[2]   = g->allNodePos[3*node3+2];

         dist = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];
         dist = 1./sqrt(dist);

         g->strandGrid[node3].normal[0] = dist*norm[0];
         g->strandGrid[node3].normal[1] = dist*norm[1];
         g->strandGrid[node3].normal[2] = dist*norm[2];

         norm[0]   = g->allNodePos[3*node4  ];
         norm[1]   = g->allNodePos[3*node4+1];
         norm[2]   = g->allNodePos[3*node4+2];

         dist = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];
         dist = 1./sqrt(dist);

         g->strandGrid[node4].normal[0] = dist*norm[0];
         g->strandGrid[node4].normal[1] = dist*norm[1];
         g->strandGrid[node4].normal[2] = dist*norm[2];

      }


   }

}
// #################################################################
//
// averageNormals
//
// Compute the average of the normals
// #################################################################
void averageNormals(GRID *g)
{

   int     i,j,temp;
   int     node1,node2,node3,node4;   
   int    *nodeCount;
   int   **indx;
   double  sumx,sumy,sumz;
   double  tang1[3],tang2[3],norm[3],dist;      
   double *tempNormal;
   
   //
   indx       = (int **)   malloc(sizeof(int *)*g->numNodePos);
   nodeCount  = (int *)    malloc(sizeof(int)*g->numNodePos);
   tempNormal = (double *) malloc(sizeof(double)*3*g->numNodePos);   

   // set all iflag and nodeCount for all nodes to 0
   for (i = 0; i < g->numNodePos; i++) nodeCount[i]=0;

   // loop through all the faces (i.e., edges)
   for (i = 0; i < g->numQuadEdge; i++)
   {
      
      // two node IDs of a given edge
      node1 = g->quadEdge[i][0];
      node2 = g->quadEdge[i][1];
      
      nodeCount[node1]++; 
      nodeCount[node2]++;      

   }

   // indx is a double pointer. So allocate memory for the
   // first level of pointer
   for (i = 0; i < g->numNodePos; i++)
      indx[i]   = (int *)    malloc(sizeof(int)*nodeCount[i]); 
      
   // set all iflag and nodeCount for all nodes to 0
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

   double facx0,facy0,facz0;
   double facx,facy,facz;
   facx0 = 1.0;
   facy0 = 1.3; //1.25;
   facz0 = 1.3; //1.25;

   // average the normals and assign to strand grids
   for (i = 0; i < g->numNodePos; i++)
   {
      if (g->flagalter[i] == 1)
      {
         facx = facx0; facy = facy0; facz = facz0;
      }
      else
      {
         facx = facy = facz = 1.;
      }


      sumx = sumy = sumz = 0.;
      for (j = 0; j < nodeCount[i]; j++)
      {
         sumx += facx*g->strandGrid[indx[i][j]].normal[0];
         sumy += facy*g->strandGrid[indx[i][j]].normal[1];
         sumz += facz*g->strandGrid[indx[i][j]].normal[2];
      }

      // average the normals
      sumx /= nodeCount[i];
      sumy /= nodeCount[i];
      sumz /= nodeCount[i];

      // normalize the new normal
      dist = sumx*sumx + sumy*sumy + sumz*sumz;
      dist = 1./sqrt(dist);

      tempNormal[3*i  ] = sumx*dist;
      tempNormal[3*i+1] = sumy*dist;
      tempNormal[3*i+2] = sumz*dist;
      
   }


   // FILE *fp = fopen("test.dat","w");
   // for (i = 0; i < g->numNodePos; i++)
   // {
   //    fprintf(fp, "%lf %lf\n",g->strandGrid[i].normal[0],tempNormal[3*i  ] );     
   //    fprintf(fp, "%lf %lf\n",g->strandGrid[i].normal[1],tempNormal[3*i+1] );     
   //    fprintf(fp, "%lf %lf\n",g->strandGrid[i].normal[2],tempNormal[3*i+2] );     
   // }
   // fclose(fp);

   // reassign the normals
   for (i = 0; i < g->numNodePos; i++)
   {
      g->strandGrid[i].normal[0] = tempNormal[3*i  ];
      g->strandGrid[i].normal[1] = tempNormal[3*i+1];
      g->strandGrid[i].normal[2] = tempNormal[3*i+2];

   }


}

// #################################################################
//
// updateVariables
//
// update outputs required for solver such as iqloops, qloops etc.
// #################################################################
void updateVariables(GRID *g)
{

   int   i,j;
   int   temp,count1;

   // allocate memoty for strandLoop and cellLoop
   g->strandLoop = (LOOP *) malloc(sizeof(LOOP)); // only strand loops
   g->cellLoop   = (LOOP *) malloc(sizeof(LOOP)); // all loops, every layer

   // length of each strand - (numStrandLayer-1)
   // total number of strands = len each strand * num quad cells
   g->strandLoop->maxLen = (g->numStrandLayer-1); // all strands equal
   g->strandLoop->totLen = (g->numStrandLayer-1)*g->numQuadConn;

   // total length of loops = ham loops per layer * total number of layers
   //                       + num strand loop 
   g->cellLoop->totLen   = g->quadLoop->totLen*(g->numStrandLayer-1)
                         + g->strandLoop->totLen;
   g->cellLoop->maxLen   = MAX(g->quadLoop->maxLen,g->strandLoop->maxLen);

   // ===============================================================
   // set the indices of the loops. Arranged as index of loops in 
   // [ layer 1, layer 2, .. , layer n-1, strands] 
   // ===============================================================   
   g->cellLoop->index     = (int *) malloc(sizeof(int)*g->cellLoop->totLen);

   //
   // accounts for the Hamiltonian paths
   //
   count1 = 0;
   temp   = pow2*g->numTriNode+1;

   for (i = 0; i < g->numStrandLayer-1; i++)
   {
      for (j = 0; j < temp; j++)
      {         
         g->cellLoop->index[count1] = g->quadLoop->index[j]
                                    + i*g->quadLoop->index[temp-1];
         count1++;
      } // j loop
   } // i loop   
   //
   // if Hybrid and strands, account for body fitted loops
   //
   if(iHybrid)
   {
      for (i = 0; i < g->numStrandLayer-1; i++)
      {
         for (j = 0; j < g->nEta-1; j++)
         {
            g->cellLoop->index[count1] = g->cellLoop->index[count1-1] 
                                       + g->nPsi;
            count1++;
         } // i loop
      } // j loop
   }
   //
   // account for the strands - indices increase by numStrandLayer   
   //
   for (i = 0; i < g->numQuadConn; i++)
   {
      g->cellLoop->index[count1] = g->cellLoop->index[count1-1] 
                                 + g->numStrandLayer;
      count1++;
   }   
   // ===============================================================
   // set the IDs of the faces of the loops. Arranged as index of 
   // loops in [ layer 1, layer 2, .. , layer n-1, strands] 
   // ===============================================================
   temp = g->quadLoop->totLen*(g->numStrandLayer-1)
        + g->numQuadConn*(g->numStrandLayer)+1;
        

   g->cellLoop->ID = (int *) malloc(sizeof(int)*temp);
   
   count1 = 0;
   //
   // account for the Hamiltonian paths at each layer
   // NOTE: Body-fitted loops are accounted for after this separately
   //
   int loopLength = g->quadLoop->totLen;
   if(iHybrid) loopLength -= g->nPsi*(g->nEta-1);

   for (i = 0; i < g->numStrandLayer-1; i++)
   {
      for (j = 0; j < loopLength; j++)
      {
         g->cellLoop->ID[count1]  = g->quadLoop->ID[j] + i*g->numQuadEdge;
         count1++;
      }
   }
   //
   // account for body-fitted loops
   //
   if(iHybrid)
   {      
      for (i = 0; i < g->numStrandLayer-1; i++)
      {
         for (j = loopLength; j < g->quadLoop->totLen; j++)
         {
            g->cellLoop->ID[count1] = g->quadLoop->ID[j] + i*g->numQuadEdge;
            count1++;
         }
      }
   } 
   //
   // account for the strands
   //
   int reverseStrands;
   reverseStrands = 1;
   if(reverseStrands)
   {
      int jumpVal,cellIDtemp;
      jumpVal                 = g->numStrandLayer-1;
      temp                    = g->numQuadEdge*(g->numStrandLayer-1);        
      cellIDtemp              = temp;

      for (i = 0; i < g->numQuadConn; i++)
      {
         g->cellLoop->ID[count1] = cellIDtemp + jumpVal;
         cellIDtemp             += g->numStrandLayer;
         count1++;
         for (j = 1; j < g->numStrandLayer; j++)
         {
            g->cellLoop->ID[count1] = g->cellLoop->ID[count1-1] - 1;
            count1++;
         }
      }
   }
   else
   {
      // account for the strands

      temp                    = g->numQuadEdge*(g->numStrandLayer-1);  
      g->cellLoop->ID[count1] = temp;   

      for (i = 0; i < g->numQuadConn; i++)
      {      
         for (j = 0; j < g->numStrandLayer; j++)
         {
            g->cellLoop->ID[count1+1] = g->cellLoop->ID[count1] + 1;
            count1++;
         }
      }
   }

   // ===============================================================
   // Create cellFaces matrix
   // ===============================================================
   g->numOctFace = g->numQuadEdge*(g->numStrandLayer-1) 
                 + g->numQuadConn*(g->numStrandLayer  );
   
   // allocate memory for octFace
   g->octFace    = (int **) malloc(sizeof(int *)*g->numOctFace);
   for (i = 0; i < g->numOctFace; i++)
      g->octFace[i] = (int *) malloc(sizeof(int)*N_CELL_FACE);

   // account for the faces formed by the Ham loops on all layers
   count1 = 0;
   for (i = 0; i < g->numStrandLayer-1; i++)
   {
      for (j = 0; j < g->numQuadEdge; j++)
      {
         // first 4 elements of octFace are node IDs
         g->octFace[count1][0] = g->quadEdge[j][0] +     i*g->numNodePos;
         g->octFace[count1][1] = g->quadEdge[j][1] +     i*g->numNodePos;
         g->octFace[count1][2] = g->quadEdge[j][1] + (i+1)*g->numNodePos;
         g->octFace[count1][3] = g->quadEdge[j][0] + (i+1)*g->numNodePos;

         // the following two are left and right Cell IDs
         g->octFace[count1][4] = g->quadEdge[j][2] +     i*g->numQuadConn;
         g->octFace[count1][5] = g->quadEdge[j][3] +     i*g->numQuadConn;

         // correction for g->octFace at edge of domain boundary
         if (g->quadEdge[j][2] == -5) g->octFace[count1][4] = -5;
         if (g->quadEdge[j][3] == -5) g->octFace[count1][5] = -5;

         // the following two are left and right face indices      
         g->octFace[count1][6] = g->quadEdge[j][4];
         g->octFace[count1][7] = g->quadEdge[j][5];



         // update counter
         count1++;

      }// j loop
   } // i loop
   //
   // take care of the strands
   //
   for (i = 0; i < g->numQuadConn; i++)
   {
      for (j = 0; j < g->numStrandLayer; j++)
      {
         // first 4 elements of octFace are node IDs
         g->octFace[count1][0] = g->quadConn[i][0] +     j*g->numNodePos;
         g->octFace[count1][1] = g->quadConn[i][1] +     j*g->numNodePos;
         g->octFace[count1][2] = g->quadConn[i][2] +     j*g->numNodePos;
         g->octFace[count1][3] = g->quadConn[i][3] +     j*g->numNodePos;

         // the following two and left and right Cell ID
         g->octFace[count1][4] = i + (j-1)*g->numQuadConn;
         g->octFace[count1][5] = i +     j*g->numQuadConn;

         // the following two are left and right face indices
         g->octFace[count1][6] = 4;
         g->octFace[count1][7] = 5;

         // Handle boundary conditions
         if (j == 0) // at surface (wall)
         {          
            g->octFace[count1][4] = -2;
            g->octFace[count1][6] = -1;
         }
         else if (j==g->numStrandLayer-1) // at free-stream
         {
            g->octFace[count1][5] = -1;
            g->octFace[count1][7] = -1;
         }

         // update counter
         count1++;

      } // j loop
   } // i loop

}

// ##################################################################
// END OF FILE
// ##################################################################
