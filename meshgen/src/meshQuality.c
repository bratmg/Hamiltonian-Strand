// ##################################################################
//
// meshQuality.c
//
// File for routines related to checking mesh quality
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "globalVariables.h"
#include "meshQuality.h"

#define NBINS  6
#define PI     3.41592653589793238462
#define R2D   57.2957795130823
#define INV90  0.01111111111111111111

void meshQuality(int gridID, GRID *g)
{
   checkQuality(gridID, g);

   outputMeshQuality(g);

}

// ##################################################################
//
// checkQuality
// Obtain the angles and aspect ratio of cells
// ##################################################################
void checkQuality(int gridID, GRID *g)
{

   int    i,j,k;
   int    n1,n2,n3,n4;
   double pts[6][3],d[5],temp;
   double theta,thMin,thMax;   
   char     filename[30],fileID[250],intStr[10]; 
   FILE   *fptr;

   // ===============================================================
   // Cell aspect ratio
   // ===============================================================

   g->meshSkewness    = (double *) malloc(sizeof(double)*g->numQuadConn);
   g->cellSizeRatio   = (double *) malloc(sizeof(double)*g->numQuadConn);
   g->numMeshSkew     = (int *)    malloc(sizeof(int)*NBINS);
   g->meshSkewnessID  = (int **)   malloc(sizeof(int *)*g->numQuadConn);

   for (i = 0; i < NBINS; i++)
   {
      g->meshSkewnessID[i] = (int *) malloc(sizeof(int)*g->numQuadConn);
      g->numMeshSkew[i]    = 0;
   }

   for (i = 0; i < g->numQuadConn; i++)
   {
      // node IDs
      n1 = g->quadConn[i][0];
      n2 = g->quadConn[i][1];
      n3 = g->quadConn[i][2];
      n4 = g->quadConn[i][3];

      // get the positions ( 4-1-2-3-4-1 )
      pts[0][0] = g->allNodePos[3*n4  ]; pts[0][1] = g->allNodePos[3*n4+1];
      pts[0][2] = g->allNodePos[3*n4+2]; 

      pts[1][0] = g->allNodePos[3*n1  ]; pts[1][1] = g->allNodePos[3*n1+1];
      pts[1][2] = g->allNodePos[3*n1+2]; 
      
      pts[2][0] = g->allNodePos[3*n2  ]; pts[2][1] = g->allNodePos[3*n2+1];
      pts[2][2] = g->allNodePos[3*n2+2]; 
      
      pts[3][0] = g->allNodePos[3*n3  ]; pts[3][1] = g->allNodePos[3*n3+1];
      pts[3][2] = g->allNodePos[3*n3+2]; 

      pts[4][0] = g->allNodePos[3*n4  ]; pts[4][1] = g->allNodePos[3*n4+1];
      pts[4][2] = g->allNodePos[3*n4+2]; 

      pts[5][0] = g->allNodePos[3*n1  ]; pts[5][1] = g->allNodePos[3*n1+1];
      pts[5][2] = g->allNodePos[3*n1+2]; 
      
      // obtain distances      
      for (j = 0; j < 5; j++)
      {
         temp = (pts[j][0] - pts[j+1][0]) * (pts[j][0] - pts[j+1][0])+
                (pts[j][1] - pts[j+1][1]) * (pts[j][1] - pts[j+1][1])+
                (pts[j][2] - pts[j+1][2]) * (pts[j][2] - pts[j+1][2]); 
                
         d[j] = sqrt(temp);         
      }

      // obtain included angles
      thMax = -1.0; thMin = 1.0e11;
      for (j = 1; j < 5; j++)
      {
         // dot product
         theta = (pts[j-1][0] - pts[j][0])*(pts[j+1][0] - pts[j][0]) +
                 (pts[j-1][1] - pts[j][1])*(pts[j+1][1] - pts[j][1]) +
                 (pts[j-1][2] - pts[j][2])*(pts[j+1][2] - pts[j][2]);  

         theta = theta/(d[j-1]*d[j]);
         theta = acos(theta)*R2D;
         
         if(theta < thMin) thMin = theta;
         if(theta > thMin) thMax = theta;
      }

      // metric for mesh skewness
      temp = INV90*MAX(thMax-90.0,90.0-thMin);

      if      (temp < 0.25) // 0-0.25 (excellent)
      {
         g->meshSkewnessID[0][g->numMeshSkew[0]] = i;
         g->numMeshSkew[0]++;          
      }
      else if (temp < 0.50) // 0.25-0.50 (good)
      {
         g->meshSkewnessID[1][g->numMeshSkew[1]] = i;
         g->numMeshSkew[1]++; 
      }
      else if (temp < 0.80) // 0.50-0.80 (acceptable)
      {
         g->meshSkewnessID[2][g->numMeshSkew[2]] = i;
         g->numMeshSkew[2]++; 
      }
      else if (temp < 0.95) // 0.80-0.95 (poor)
      {
         g->meshSkewnessID[3][g->numMeshSkew[3]] = i;
         g->numMeshSkew[3]++; 
      }
      else if (temp < 0.99) // 0.95-0.99 (sliver)
      {
         g->meshSkewnessID[4][g->numMeshSkew[4]] = i;
         g->numMeshSkew[4]++; 
      }
      else                  // 0.99-1.00 (degenerate)
      {         
         g->meshSkewnessID[5][g->numMeshSkew[5]] = i;
         g->numMeshSkew[5]++;
      }

      g->meshSkewness[i] = temp;


   } // i loop


   // write out cell connectivity information to file
   strcpy(fileID,folderName);
   strcat(fileID,"cellQuality.tecdat");
   fptr = fopen(fileID,"w");
   for (i = 0; i < NBINS; i++)
   {
      if(g->numMeshSkew[i]>0)
      {
         fprintf(fptr,"TITLE = \"BIN ID %d \" \n",i);
         fprintf(fptr,"VARIABLES=\"X\",\"Y\",\"Z\" \n");
         fprintf(fptr,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",
            g->numNodePos,g->numMeshSkew[i]);

         // write the position information
         k = 0;
         for (j = 0; j <  g->numNodePos; j++)
         {
            fprintf(fptr, "%lf %lf %lf \n", 
               g->allNodePos[k],g->allNodePos[k+1],g->allNodePos[k+2]);
            k += 3;
         }

         // write the connectivity information   
         for (j = 0; j <  g->numMeshSkew[i]; j++)
         {
            k = g->meshSkewnessID[i][j];
            fprintf(fptr, "%d %d %d %d \n", 
               g->quadConn[k][0]+1,g->quadConn[k][1]+1,
               g->quadConn[k][2]+1,g->quadConn[k][3]+1);
         }
      }
   }

   fclose(fptr);

   // write cell quality information for histogram
   strcpy(fileID,folderName);
   strcat(fileID,"histogram.dat");
   fptr = fopen(fileID,"w");
   for (i = 0; i < g->numQuadConn; i++)
      fprintf(fptr,"%lf \n",g->meshSkewness[i]);

   fclose(fptr);


   // ===============================================================
   // Cell size quality (as a fraction of the original triangle from
   // which it was formed)
   // ===============================================================
   g->triArea  = (double *) malloc(sizeof(double)*g->numTriangle);
   g->quadArea = (double *) malloc(sizeof(double)*g->numQuadConn);

   double a,b,c,e,s,p,q,area;

   // compute area of all triangles
   for (i = 0; i < g->numTriangle; i++)
   {
      n1 = g->triConn[3*i    ];
      n2 = g->triConn[3*i + 1];
      n3 = g->triConn[3*i + 2];

      a = (g->allNodePos[3*n1  ] - g->allNodePos[3*n2  ])*
          (g->allNodePos[3*n1  ] - g->allNodePos[3*n2  ]) + 
          (g->allNodePos[3*n1+1] - g->allNodePos[3*n2+1])*
          (g->allNodePos[3*n1+1] - g->allNodePos[3*n2+1]) + 
          (g->allNodePos[3*n1+2] - g->allNodePos[3*n2+2])*
          (g->allNodePos[3*n1+2] - g->allNodePos[3*n2+2]);
      a = sqrt(a);

      b = (g->allNodePos[3*n1  ] - g->allNodePos[3*n3  ])*
          (g->allNodePos[3*n1  ] - g->allNodePos[3*n3  ]) + 
          (g->allNodePos[3*n1+1] - g->allNodePos[3*n3+1])*
          (g->allNodePos[3*n1+1] - g->allNodePos[3*n3+1]) + 
          (g->allNodePos[3*n1+2] - g->allNodePos[3*n3+2])*
          (g->allNodePos[3*n1+2] - g->allNodePos[3*n3+2]);
      b = sqrt(b);

      c = (g->allNodePos[3*n3  ] - g->allNodePos[3*n2  ])*
          (g->allNodePos[3*n3  ] - g->allNodePos[3*n2  ]) + 
          (g->allNodePos[3*n3+1] - g->allNodePos[3*n2+1])*
          (g->allNodePos[3*n3+1] - g->allNodePos[3*n2+1]) + 
          (g->allNodePos[3*n3+2] - g->allNodePos[3*n2+2])*
          (g->allNodePos[3*n3+2] - g->allNodePos[3*n2+2]);
      c = sqrt(c);

      // semiperimeter
      s = 0.5*(a+b+c);

      g->triArea[i] = sqrt(s*(s-a)*(s-b)*(s-c));

   } // i loop

   // compute area of quadrilaterals
   for (i = 0; i < g->numQuadConn; i++)
   {
      n1 = g->quadConn[i][0];
      n2 = g->quadConn[i][1];
      n3 = g->quadConn[i][2];
      n4 = g->quadConn[i][3];

      a = (g->allNodePos[3*n1  ] - g->allNodePos[3*n2  ])*
          (g->allNodePos[3*n1  ] - g->allNodePos[3*n2  ]) + 
          (g->allNodePos[3*n1+1] - g->allNodePos[3*n2+1])*
          (g->allNodePos[3*n1+1] - g->allNodePos[3*n2+1]) + 
          (g->allNodePos[3*n1+2] - g->allNodePos[3*n2+2])*
          (g->allNodePos[3*n1+2] - g->allNodePos[3*n2+2]);
      a = sqrt(a);

      b = (g->allNodePos[3*n2  ] - g->allNodePos[3*n3  ])*
          (g->allNodePos[3*n2  ] - g->allNodePos[3*n3  ]) + 
          (g->allNodePos[3*n2+1] - g->allNodePos[3*n3+1])*
          (g->allNodePos[3*n2+1] - g->allNodePos[3*n3+1]) + 
          (g->allNodePos[3*n2+2] - g->allNodePos[3*n3+2])*
          (g->allNodePos[3*n2+2] - g->allNodePos[3*n3+2]);
      b = sqrt(b);

      c = (g->allNodePos[3*n3  ] - g->allNodePos[3*n4  ])*
          (g->allNodePos[3*n3  ] - g->allNodePos[3*n4  ]) + 
          (g->allNodePos[3*n3+1] - g->allNodePos[3*n4+1])*
          (g->allNodePos[3*n3+1] - g->allNodePos[3*n4+1]) + 
          (g->allNodePos[3*n3+2] - g->allNodePos[3*n4+2])*
          (g->allNodePos[3*n3+2] - g->allNodePos[3*n4+2]);
      c = sqrt(c);

      e = (g->allNodePos[3*n4  ] - g->allNodePos[3*n1  ])*
          (g->allNodePos[3*n4  ] - g->allNodePos[3*n1  ]) + 
          (g->allNodePos[3*n4+1] - g->allNodePos[3*n1+1])*
          (g->allNodePos[3*n4+1] - g->allNodePos[3*n1+1]) + 
          (g->allNodePos[3*n4+2] - g->allNodePos[3*n1+2])*
          (g->allNodePos[3*n4+2] - g->allNodePos[3*n1+2]);
      e = sqrt(e);

      p = (g->allNodePos[3*n1  ] - g->allNodePos[3*n3  ])*
          (g->allNodePos[3*n1  ] - g->allNodePos[3*n3  ]) + 
          (g->allNodePos[3*n1+1] - g->allNodePos[3*n3+1])*
          (g->allNodePos[3*n1+1] - g->allNodePos[3*n3+1]) + 
          (g->allNodePos[3*n1+2] - g->allNodePos[3*n3+2])*
          (g->allNodePos[3*n1+2] - g->allNodePos[3*n3+2]);
      p = sqrt(p);

      // q = (g->allNodePos[3*n2  ] - g->allNodePos[3*n4  ])*
      //     (g->allNodePos[3*n2  ] - g->allNodePos[3*n4  ]) + 
      //     (g->allNodePos[3*n2+1] - g->allNodePos[3*n4+1])*
      //     (g->allNodePos[3*n2+1] - g->allNodePos[3*n4+1]) + 
      //     (g->allNodePos[3*n2+2] - g->allNodePos[3*n4+2])*
      //     (g->allNodePos[3*n2+2] - g->allNodePos[3*n4+2]);
      // q = sqrt(q);      


      s = 0.5*(a+b+p);

      g->quadArea[i] = sqrt(s*(s-a)*(s-b)*(s-p));

      s = 0.5*(c+e+p);

      g->quadArea[i] += sqrt(s*(s-c)*(s-e)*(s-p));


      g->cellSizeRatio[i] = 3.0*pow4*g->quadArea[i]
                          /(g->triArea[g->quad2triList[i]]);

   } // i loop

   // write cell quality information for histogram
   strcpy(fileID,folderName);
   strcat(fileID,"cellsize.dat");
   fptr = fopen(fileID,"w");
   for (i = 0; i < g->numQuadConn; i++)
      fprintf(fptr,"%lf \n",g->cellSizeRatio[i]);

   fclose(fptr);



}

// ##################################################################
//
// outputMeshQuality
//
// Outputs the statistics to a file.
//
// My only request to those who add features in the future is to 
// maintain the asthetic sense of the file, so as to ensure readability
// ##################################################################
void outputMeshQuality(GRID *g)
{

   int   i;
   char  fileID[250];
   strcpy(fileID,folderName);
   strcat(fileID,"statistics.dat");
   FILE *fptr = fopen(fileID,"w");

   fprintf(fptr,"#####################################################################\n");
   fprintf(fptr,"#                                                                   #\n");
   fprintf(fptr,"#   STATISTICS: MESH GEN WITH HAMILTONIAN LOOPS AND STRAND GRIDS    #\n");
   fprintf(fptr,"#                                                                   #\n");
   fprintf(fptr,"#####################################################################\n");  
   
   // ===============================================================
   fprintf(fptr,"\n General statistics:\n");
   // ===============================================================
   fprintf(fptr,"   Number of input nodes:        %d\n",g->numTriNode);
   fprintf(fptr,"   Number total nodes per layer: %d\n",g->numNodePos);   

   if(iStrand) fprintf(fptr,"   Number of strand layers:      %d\n",g->numStrandLayer);
   fprintf(fptr,"   Number of quads:              %d\n",g->numQuadConn);

   if(iStrand) fprintf(fptr,"   Number of cells:              %d\n",
                       g->numQuadConn*(g->numStrandLayer-1));
   
   fprintf(fptr,"   Number of quad edges:         %d\n",g->numQuadEdge);
   if(iStrand) fprintf(fptr,"   Number of oct faces:          %d\n",g->numOctFace);

	fprintf(fptr,"   Number of Ham loops:          %d\n",pow2*g->numTriNode+1);
	fprintf(fptr,"   Longest Ham loop:             %d\n",g->quadLoop->maxLen);
   fprintf(fptr,"\n   Total number of colours:      %d\n",g->maxCol+1);
   for (i = 0; i < g->maxCol+1; i++)
      fprintf(fptr,"    - Loops of colour %d:         %d\n",i+1,g->numCol[i]);
   // ===============================================================
   fprintf(fptr,"\n 2d mesh quality statistics:\n");
   // ===============================================================
   fprintf(fptr,"   Cell skewness (0.00-0.25):   %d\n",g->numMeshSkew[0]);
   fprintf(fptr,"   Cell skewness (0.25-0.50):   %d\n",g->numMeshSkew[1]);
   fprintf(fptr,"   Cell skewness (0.50-0.80):   %d\n",g->numMeshSkew[2]);
   fprintf(fptr,"   Cell skewness (0.80-0.95):   %d\n",g->numMeshSkew[3]);
   fprintf(fptr,"   Cell skewness (0.95-0.99):   %d\n",g->numMeshSkew[4]);
   fprintf(fptr,"   Cell skewness (0.99-1.00):   %d\n",g->numMeshSkew[5]);

   if(iStrand)
   {
      // ============================================================
      //printf(fptr,"\n 3d mesh quality statistics:\n");
      // ============================================================



   }

   // ===============================================================
   //fprintf(fptr,"\n Memory statistics:\n");
   // ===============================================================

   fprintf(fptr,"\n");
   fprintf(fptr,"#####################################################################\n");
   fprintf(fptr,"#                                                                   #\n");
   fprintf(fptr,"#                          END OF FILE                              #\n");
   fprintf(fptr,"#                                                                   #\n");
   fprintf(fptr,"#####################################################################\n");

   fclose(fptr);

}


// ##################################################################
// END OF FILE
// ##################################################################
