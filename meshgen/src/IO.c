// ##################################################################
//
// IO.c
//
// File that contains the input/output routines
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "globalVariables.h"
#include "IO.h"
#include "meshtype.h"

// ##################################################################
// WRITE AN ARRAY TO FILE - ONLY FOR TESTING
// ##################################################################
void writearrayINT(const int *array, const int n)
{
   int i;
   
   FILE * fp = fopen("test.dat","w");   

   for(i=0; i<n; i++)   
      fprintf(fp,"%d %d\n",i,array[i]);
   
   fclose(fp);
}


void writearrayDUB(const double *array, const int n)
{
   int i;
   
   FILE * fp = fopen("test.dat","w");

   for(i=0; i<n; i++)   
      fprintf(fp,"%d %lf\n",i,array[i]);


   fclose(fp);
}

void writemultiarrayINT(const int **array, const int m, const int n)
{
   int i,j;
   
   FILE * fp = fopen("test.dat","w");

   for(i=0; i<m; i++)   
      for (j=0; j<n; j++)
         fprintf(fp,"%d %d %d\n",i,j,array[i][j]);


   fclose(fp);
}

void writemultiarrayDUB(const double **array, const int m, const int n)
{
   int i,j;
   
   FILE * fp = fopen("test.dat","w");

   for(i=0; i<m; i++)   
      for (j=0; j<n; j++)
         fprintf(fp,"%d %d %lf\n",i,j,array[i][j]);


   fclose(fp);
}

// ##################################################################
//
// welcome and thanks
//
// ##################################################################
void welcome()
{
   printf("#####################################################################\n");
   printf("#                                                                   #\n");
   printf("#     MESH GENERATION USING HAMILTONIAN PATHS AND STRAND GRIDS      #\n");
   printf("#                                                                   #\n");
   printf("#####################################################################\n\n");

}

void thanks()
{
   printf("\n#####################################################################\n");
   printf("#                                                                   #\n");
   printf("#                             FIN.                                  #\n");
   printf("#                                                                   #\n");
   printf("#####################################################################\n");

}

// ##################################################################
// readInputs 
//  reads from input.meshgen and from the node and connectivity data
// ##################################################################
void readInputs()
{
	FILE * fptr;
   
   int    i,j,k;
   int    nLayer,qLevel;
   char   c,buffer[100];
   char   nodeData[250];
   char   connData[250];
   double growth,initLen;

   // ===============================================================
   // Read input.meshgen
   // ===============================================================
   if( (fptr = fopen("input.meshgen","r")) == NULL)
   {
      printf("Error: Input file 'input.meshgen' missing.\n");
      exit(1);
   }

   // skip garbage lines
   for (i = 0 ; i < 5; i++)
      while((c=fgetc(fptr))!='\n'); 

   fscanf(fptr,"num grid      = %d\n",&nGrid);
   fscanf(fptr,"node file     = %s\n",nodeData);
   fscanf(fptr,"conn file     = %s\n",connData);
   fscanf(fptr,"flush surface = %s\n",surfaceType);
   fscanf(fptr,"quad level    = %d\n",&qLevel);
   fscanf(fptr,"num smooth    = %d\n",&numSmooth);
   fscanf(fptr,"if strands    = %d\n",&iStrand);
   fscanf(fptr,"strand layers = %d\n",&nLayer);
   fscanf(fptr,"init spacing  = %lf\n",&initLen);
   fscanf(fptr,"mesh growth   = %lf\n",&growth);
   fclose(fptr);

   if(iStrand && nLayer < 2) nLayer = 2;

   // initialize the number of grids
   g = (GRID *) malloc(sizeof(GRID)*nGrid);


   // ===============================================================
   // Loop through number of grids
   // ===============================================================
   for (i = 0; i < nGrid; i++) // loop through total number of grids
   {
      // ============================================================
      // Read coordinate data
      // ============================================================
      if( (fptr = fopen(nodeData,"r")) == NULL)
      {
         printf("Error: Node data file missing. \n");
         exit(1);
      }

      // obtain the number of node
      fscanf(fptr, "%d\n", &g[i].numTriNode);

      // allocate space in nodePosTri
      g[i].nodePosTri = (double *) malloc(sizeof(double)*3*g[i].numTriNode);


      k = 0;
      for (j = 0; j < g[i].numTriNode; j++)
      {
         fscanf(fptr, "%lf %lf %lf \n",
            &g[i].nodePosTri[k],&g[i].nodePosTri[k+1],&g[i].nodePosTri[k+2]);
      
         k += 3;
      }   

      fclose(fptr); // close nodeData


      // ============================================================ 
      // Read connectivity data
      // ============================================================
      if( (fptr = fopen(connData,"r")) == NULL)
      {
         printf("Error: Connectivity data file missing. \n");
         exit(1);
      }


      // obtain the number of node
      fscanf(fptr, "%d\n", &g[i].numTriangle);

      // allocate space in nodePosTri
      g[i].triConn = (int *) malloc(sizeof(int)*3*g[i].numTriangle);

      k = 0;
      for (j = 0; j < g[i].numTriangle; j++)
      {
         fscanf(fptr, "%d %d %d \n",
            &g[i].triConn[k],&g[i].triConn[k+1],&g[i].triConn[k+2]);


         // Because the connectivity is created using a 1 based indexing
         // language (Matlab), subtract all by 1 to make it amenable to C
         g[i].triConn[k  ]--;
         g[i].triConn[k+1]--;
         g[i].triConn[k+2]--;
         k += 3;
      }   

      fclose(fptr); // close nodeData

      // ============================================================
      // Set variables related to strand meshes
      // ============================================================
      g[i].numStrandLayer = nLayer;
      g[i].initMeshLen    = initLen;
      g[i].meshGrowth     = growth;
      
      // quad division level
      if(qLevel < 1)
      {
         printf("#meshgen: quad level of %d is incompatible.\n",qLevel);
         printf("#meshgen: Setting quad level to 1");
         qLevel = 1;
      }
      g[i].quadLevel    = qLevel;
      g[i].nVertPerSide = pow(2,g->quadLevel+1)-1;
      g[i].nEdgePerSide = g[i].nVertPerSide+1;
      g[i].nIntPts      = 3*pow(2,2*g->quadLevel-2) +
                          3*pow(2,g->quadLevel-1) + 1;


		// set global variables pow2 and pow4
		pow2 = pow(2,g[i].quadLevel);
		pow4 = pow(4,g[i].quadLevel);

   } // end i for number of grid


}

// ##################################################################
//
// writeTecplot
//  reads from input.meshgen and from the node and connectivity data
//
// ##################################################################
void writeTecplot(GRID *g)
{

   int      i,j,k,kk,temp;
   int      ii,i1,i2,jj;
   int      ie,iA,iB,iC,iD;
   int      count1,count2,lenCount;
   int     *edgeConn = (int *) malloc(sizeof(int)*2*g->numQuadEdge);
   char     filename[30]; 
   double   dist;
   double  *edgePts  = (double *)  malloc(sizeof(double)*3*g->quadLoop->totLen);   
   double **lenLoop  = (double **) malloc(sizeof(double *)*pow2*g->numTriNode);
   FILE    *fptr;

   // ===============================================================
   // TRIANGLES.DAT
   // ===============================================================   
   fptr = fopen("./output/triangles.tecdat","w");
   fprintf(fptr,"VARIABLES=\"X\",\"Y\",\"Z\" \n");
   fprintf(fptr,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",
      g->numTriNode,g->numTriangle);

   // write the position information
   k = 0;
   for (j = 0; j <  g->numTriNode; j++)
   {
      fprintf(fptr, "%lf %lf %lf \n", 
         g->nodePosTri[k],g->nodePosTri[k+1],g->nodePosTri[k+2]);
      k += 3;
   }

   // write the connectivity information
   k = 0;
   for (j = 0; j <  g->numTriangle; j++)
   {
      fprintf(fptr, "%d %d %d \n", 
         g->triConn[k]+1,g->triConn[k+1]+1,g->triConn[k+2]+1);
      k += 3;
   }

   fclose(fptr);

   // ===============================================================
   // QUADS.DAT
   // ===============================================================   
   fptr = fopen("./output/quads.tecdat","w");
   fprintf(fptr,"VARIABLES=\"X\",\"Y\",\"Z\" \n");
   fprintf(fptr,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",
      g->numNodePos,g->numQuadConn);

   // write the position information
   k = 0;
   for (j = 0; j <  g->numNodePos; j++)
   {
      fprintf(fptr, "%lf %lf %lf \n", 
         g->allNodePos[k],g->allNodePos[k+1],g->allNodePos[k+2]);
      k += 3;
   }

   // write the connectivity information   
   for (j = 0; j <  g->numQuadConn; j++)
   {
      fprintf(fptr, "%d %d %d %d \n", 
         g->quadConn[j][0]+1,g->quadConn[j][1]+1,
         g->quadConn[j][2]+1,g->quadConn[j][3]+1);
   }
   fclose(fptr);

   // ===============================================================
   // COORD.DAT
   // ===============================================================   
   fptr = fopen("./output/coord.dat","w");
   if(iStrand) 
   {
      // create node positions

      count1 = g->numNodePos*g->numStrandLayer;

   }
   else
   {
      count1 = g->numNodePos;
   }

   fprintf(fptr,"%d\n",count1);
   
   // write the position information   
   // I believe this write is pretty inefficient as the
   // indexing is upside down. Changes the pointer faster than
   // the pos index values.
   if(iStrand)
   {   
      
      g->cellNodePos = (double *) malloc(sizeof(double)*3*count1);

      k = 0;
      for (i = 0; i < g->numStrandLayer; i++)
      {
         for (j = 0; j < g->numNodePos; j++)
         {
            g->cellNodePos[k  ] = g->strandGrid[j].pos[3*i  ];
            g->cellNodePos[k+1] = g->strandGrid[j].pos[3*i+1];
            g->cellNodePos[k+2] = g->strandGrid[j].pos[3*i+2];            

            fprintf(fptr,"%lf %lf %lf \n",
               g->cellNodePos[k  ],g->cellNodePos[k+1],g->cellNodePos[k+2]);

            k+=3;
         } // j loop
      } // i loop   

   }
   else
   {
      k = 0;
      for (j = 0; j <  g->numNodePos; j++)
      {
         fprintf(fptr, "%lf %lf %lf \n", 
            g->allNodePos[k],g->allNodePos[k+1],g->allNodePos[k+2]);
         k += 3;
      }      
   }
   fclose(fptr);

   // ===============================================================
   // CONN.DAT
   // ===============================================================   
   if(iStrand)
   {
      fptr = fopen("./output/conn.dat","w");
      fprintf(fptr,"%d\n",g->numQuadConn*(g->numStrandLayer-1));   
      for (i = 0; i < g->numStrandLayer-1; i++)
      {
         for (j = 0; j <  g->numQuadConn; j++)
         {
            fprintf(fptr, "%d %d %d %d %d %d %d %d\n", 
               i*g->numNodePos+(g->quadConn[j][0]+1),
               i*g->numNodePos+(g->quadConn[j][1]+1),
               i*g->numNodePos+(g->quadConn[j][2]+1),
               i*g->numNodePos+(g->quadConn[j][3]+1),
           (i+1)*g->numNodePos+(g->quadConn[j][0]+1),
           (i+1)*g->numNodePos+(g->quadConn[j][1]+1),
           (i+1)*g->numNodePos+(g->quadConn[j][2]+1),
           (i+1)*g->numNodePos+(g->quadConn[j][3]+1));
         }

      }
      fclose(fptr);
   }
   else
   {
      fptr = fopen("./output/conn.dat","w");
      fprintf(fptr,"%d\n",g->numQuadConn);
      
      // write the connectivity information   
      for (j = 0; j <  g->numQuadConn; j++)
      {
         fprintf(fptr, "%d %d %d %d \n", 
            g->quadConn[j][0],g->quadConn[j][1],
            g->quadConn[j][2],g->quadConn[j][3]);
      }
      fclose(fptr);
   }

   // ===============================================================
   // QEDGES.DAT / OFACES.DAT
   // ===============================================================   
   if (iStrand)
   {
      fptr = fopen("./output/ofaces.dat","w");
      fprintf(fptr,"%d\n",g->numOctFace);
            
      for (j = 0; j <  g->numOctFace; j++)
      {
         fprintf(fptr, "%d %d %d %d %d %d %d %d\n", 
            g->octFace[j][0],g->octFace[j][1],g->octFace[j][2],
            g->octFace[j][3],g->octFace[j][4],g->octFace[j][5],
            g->octFace[j][6],g->octFace[j][7]);
      }
      fclose(fptr);
   }
   else
   {
      fptr = fopen("./output/qedges.dat","w");
      fprintf(fptr,"%d\n",g->numQuadEdge);

      for (j = 0; j <  g->numQuadEdge; j++)
      {
         fprintf(fptr, "%d %d %d %d %d %d\n", 
            g->quadEdge[j][0],g->quadEdge[j][1],g->quadEdge[j][2],
            g->quadEdge[j][3],g->quadEdge[j][4],g->quadEdge[j][5]);
      }
      fclose(fptr);  
   }

   // ===============================================================
   // QLOOPS.DAT
   // ===============================================================   
   fptr = fopen("./output/qloops.dat","w");
   if(iStrand)
   {
      // total number of faces =   faces in each layer * (numLayer-1)
      //                    + number of strands * number of faces in each strand 
      temp = g->quadLoop->totLen*(g->numStrandLayer-1)
                          + g->numQuadConn*(g->numStrandLayer);
      fprintf(fptr,"%d\n", temp);

      for (j = 0; j < temp; j++)
         fprintf(fptr, "%d\n", g->cellLoop->ID[j]);
   }
   else
   {
      fprintf(fptr,"%d\n",g->quadLoop->totLen);
         
      for (j = 0; j < g->quadLoop->totLen; j++)
         fprintf(fptr, "%d\n", g->quadLoop->ID[j]);
   }

   fclose(fptr);

   // ===============================================================
   // IQLOOPS.DAT
   // ===============================================================      
   fptr = fopen("./output/iqloops.dat","w");
   if(iStrand) 
   {
      // understand the reason for +2
      temp = (pow2*g->numTriNode+1)*(g->numStrandLayer-1)+g->numQuadConn;
      fprintf(fptr,"%d\n",temp);
      for (j = 0; j < temp; j++)
         fprintf(fptr, "%d\n", g->cellLoop->index[j]);
   }
   else
   {
      temp = pow2*g->numTriNode+1;
      fprintf(fptr,"%d\n",temp);
         
      for (j = 0; j < temp; j++)
         fprintf(fptr, "%d\n", g->quadLoop->index[j]);
         
   }
   
   fclose(fptr);

   // ===============================================================
   // NCOLORS.DAT
   // ===============================================================

   fptr = fopen("./output/ncolors.dat","w");
   if(iStrand)
      fprintf(fptr,"%d\n",g->maxCol+2);
   else
      fprintf(fptr,"%d\n",g->maxCol+1);
   
   // write the connectivity information   
   for (j = 0; j <= g->maxCol; j++)
      fprintf(fptr, "%d\n", g->numCol[j]);
   
   if (iStrand) fprintf(fptr, "%d\n", g->numQuadConn);


   fclose(fptr);

   // ===============================================================
   // TECPLOT LOOPS.DAT/ STRANDS_LOOPS.DAT
   // ===============================================================
   if(iStrand)
   {
      int ktemp;
      temp  = pow2*g->numQuadEdge*(g->numStrandLayer-1)
                + g->numStrandLayer*g->numQuadConn + 1;
      int    *faceConn = (int *)     malloc(sizeof(int)*temp);
      double *facePts  = (double *)  malloc(sizeof(double)*3*g->cellLoop->totLen);   

      k     = 0;
      ktemp = 0;
      for (i = 0; i <= g->maxCol; i++)      
      {

         printf("#meshgen: Preparing Tecplot file for loops %d ...\n",i+1);
         sprintf(filename,"./output/loops%d.tecdat", i+1);      
         fptr = fopen(filename,"w");

         count1 = 0;
         count2 = 0;
         kk     = 0;
         temp   = pow2*g->numCol[i];  

         for (jj = 0; jj < g->numStrandLayer-1; jj++)         
         {
            
            k = ktemp;            
            for (j = 0; j < temp; j++)
            {         
               i1 = g->cellLoop->index[k   + jj*(pow2*g->numTriNode+1)];
               i2 = g->cellLoop->index[k+1 + jj*(pow2*g->numTriNode+1)];                                 
               k++;

               for (ii = i1; ii < i2; ii++)
               {  
                  ie = g->cellLoop->ID[ii];
                  iA = g->octFace[ie][0];
                  iB = g->octFace[ie][1];
                  iC = g->octFace[ie][2];
                  iD = g->octFace[ie][3];

                  facePts[3*count1  ] = 0.25*(g->cellNodePos[3*iA  ] + g->cellNodePos[3*iB  ] 
                                      +       g->cellNodePos[3*iC  ] + g->cellNodePos[3*iD  ]);
                  facePts[3*count1+1] = 0.25*(g->cellNodePos[3*iA+1] + g->cellNodePos[3*iB+1]
                                      +       g->cellNodePos[3*iC+1] + g->cellNodePos[3*iD+1]);
                  facePts[3*count1+2] = 0.25*(g->cellNodePos[3*iA+2] + g->cellNodePos[3*iB+2]
                                      +       g->cellNodePos[3*iC+2] + g->cellNodePos[3*iD+2]);

                  count1++;
                  if (ii < i2-1)
                  {               
                     faceConn[2*count2  ] = kk;
                     faceConn[2*count2+1] = kk+1;
                     count2++;
                  }
                  kk++;
               } // ii loop
            } // j loop
         } // jj loop   
         ktemp += temp;  

         fprintf(fptr,"VARIABLES=\"X\",\"Y\",\"Z\" \n");
         fprintf(fptr,"ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG\n",
         count1,count2);

         for (j = 0; j < count1; j++)
            fprintf(fptr,"%lf %lf %lf \n",facePts[3*j],facePts[3*j+1],facePts[3*j+2]);

         for (j = 0; j < count2; j++)
            fprintf(fptr,"%d %d \n",faceConn[2*j]+1,faceConn[2*j+1]+1);

         fclose(fptr);
         
      }

      // NSTRANDS.DAT
      fptr   = fopen("./output/nstrands.dat","w");
      fprintf(fptr,"%d\n",g->numStrandLayer);
      fclose(fptr);

      // STRANDS.DAT
      printf("#meshgen: Preparing Tecplot file for strands ...\n");         
      fptr   = fopen("./output/strands.tecdat","w");
      
      k      = (pow2*g->numTriNode+1)*(g->numStrandLayer-1)-1;

      count1 = 0;
      count2 = 0;
      kk     = 0;
      temp   = g->numQuadConn;

      for (j = 0; j < temp; j++)
      {         
         i1 = g->cellLoop->index[k  ];
         i2 = g->cellLoop->index[k+1];         
         k++;
         
         for (ii = i1; ii < i2; ii++)
         {  
            ie = g->cellLoop->ID[ii];            
            iA = g->octFace[ie][0];
            iB = g->octFace[ie][1];
            iC = g->octFace[ie][2];
            iD = g->octFace[ie][3];

            facePts[3*count1  ] = 0.25*(g->cellNodePos[3*iA  ] + g->cellNodePos[3*iB  ] 
                                +       g->cellNodePos[3*iC  ] + g->cellNodePos[3*iD  ]);
            facePts[3*count1+1] = 0.25*(g->cellNodePos[3*iA+1] + g->cellNodePos[3*iB+1]
                                +       g->cellNodePos[3*iC+1] + g->cellNodePos[3*iD+1]);
            facePts[3*count1+2] = 0.25*(g->cellNodePos[3*iA+2] + g->cellNodePos[3*iB+2]
                                +       g->cellNodePos[3*iC+2] + g->cellNodePos[3*iD+2]);

            count1++;
            if (ii < i2-1)
            {               
               faceConn[2*count2  ] = kk;
               faceConn[2*count2+1] = kk+1;
               count2++;
            }
            kk++;
         } // ii loop
      } // j loop
      

      fprintf(fptr,"VARIABLES=\"X\",\"Y\",\"Z\" \n");
      fprintf(fptr,"ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG\n",
      count1,count2);

      for (j = 0; j < count1; j++)
         fprintf(fptr,"%lf %lf %lf \n",facePts[3*j],facePts[3*j+1],facePts[3*j+2]);

      for (j = 0; j < count2; j++)
         fprintf(fptr,"%d %d \n",faceConn[2*j]+1,faceConn[2*j+1]+1);

      fclose(fptr);
       

   }
   else
   {
     for (i = 0; i < pow2*g->numTriNode; i++)
         lenLoop[i] = (double *) malloc(sizeof(double)*g->quadLoop->maxLen);

      for (i = 0; i < pow2*g->numTriNode; i++)
         for (j = 0; j < g->quadLoop->maxLen; j++)
            lenLoop[i][j] = -1.0;

      k = 0;
      for (i = 0; i <= g->maxCol; i++)
      {

         printf("#meshgen: Preparing Tecplot file for loops %d ...\n",i+1);
         sprintf(filename,"./output/loops%d.tecdat", i+1);      
         fptr = fopen(filename,"w");

    
         count1 = 0;
         count2 = 0;
         kk     = 0;
         temp   = pow2*g->numCol[i];  
         
         for (j = 0; j < temp; j++)
         {         
            i1 = g->quadLoop->index[k];
            i2 = g->quadLoop->index[k+1];
            k++;         

            lenCount = 0; 
            for (ii = i1; ii < i2; ii++)
            {  
               ie = g->quadLoop->ID[ii];
               iA = g->quadEdge[ie][0];
               iB = g->quadEdge[ie][1];

               edgePts[3*count1  ] = 0.5*(g->allNodePos[3*iA  ] + g->allNodePos[3*iB  ]);
               edgePts[3*count1+1] = 0.5*(g->allNodePos[3*iA+1] + g->allNodePos[3*iB+1]);
               edgePts[3*count1+2] = 0.5*(g->allNodePos[3*iA+2] + g->allNodePos[3*iB+2]);


               // // compute L2 distance
               // if (ii != i1)
               // {
               //    dist = (edgePts[3*count1  ]-edgePts[3*(count1-1)  ])*
               //           (edgePts[3*count1  ]-edgePts[3*(count1-1)  ]) +
               //           (edgePts[3*count1+1]-edgePts[3*(count1-1)+1])*
               //           (edgePts[3*count1+1]-edgePts[3*(count1-1)+1]) +
               //           (edgePts[3*count1+2]-edgePts[3*(count1-1)+2])*
               //           (edgePts[3*count1+2]-edgePts[3*(count1-1)+2]);

               //    lenLoop[k-1][lenCount]=sqrt(dist);
               //    lenCount++;
               // }

               count1++;                 
               if (ii < i2-1)
               {               
                  edgeConn[2*count2  ] = kk;
                  edgeConn[2*count2+1] = kk+1;
                  count2++;
               }
               kk++;
            } // ii loop
         } // j loop      

         fprintf(fptr,"VARIABLES=\"X\",\"Y\",\"Z\" \n");
         fprintf(fptr,"ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG\n",
         count1,count2);

         for (j = 0; j < count1; j++)
            fprintf(fptr,"%lf %lf %lf \n",edgePts[3*j],edgePts[3*j+1],edgePts[3*j+2]);

         for (j = 0; j < count2; j++)
            fprintf(fptr,"%d %d \n",edgeConn[2*j]+1,edgeConn[2*j+1]+1);

         fclose(fptr);

         // // write the length of the loops
         // sprintf(filename,"./output/lenLoops%d.tecdat", i+1);      
         // fptr = fopen(filename,"w");
         // fprintf(fptr,"%d\n",2*g->numCol[i]);
         // fprintf(fptr,"%d\n",g->quadLoop->maxLen);
         // for (ii = 0; ii < 2*g->numTriNode; ii++)
         //    for (j = 0; j < g->quadLoop->maxLen; j++)
         //       fprintf(fptr,"%lf\n",lenLoop[ii][j]);

         // fclose(fptr);
      }

   }

   // ===============================================================
   // BRICK.DAT (only if strands are required)
   // ===============================================================
   if(iStrand)
   {      
      sprintf(filename,"./output/brick.tecdat");
      fptr = fopen(filename,"w");

      count1 = g->numStrandLayer*g->numNodePos;      // num nodes
      count2 = g->numQuadConn*(g->numStrandLayer-1); // num edges

      fprintf(fptr,"VARIABLES=\"X\",\"Y\",\"Z\" \n");
      fprintf(fptr,"ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FEBRICK\n",
      count1,count2);
      
      for (i = 0; i < g->numStrandLayer; i++)
      {
         for (j = 0; j < g->numNodePos; j++)
         {
            fprintf(fptr,"%lf %lf %lf \n",
               g->strandGrid[j].pos[3*i  ],
               g->strandGrid[j].pos[3*i+1],
               g->strandGrid[j].pos[3*i+2]);
         } // j loop
      } // i loop      

      // connectivity information
      for (i = 0; i < g->numStrandLayer-1; i++)
      {
         for (j = 0; j <  g->numQuadConn; j++)
         {
            fprintf(fptr, "%d %d %d %d %d %d %d %d\n", 
               i*g->numNodePos+(g->quadConn[j][0]+1),
               i*g->numNodePos+(g->quadConn[j][1]+1),
               i*g->numNodePos+(g->quadConn[j][2]+1),
               i*g->numNodePos+(g->quadConn[j][3]+1),
           (i+1)*g->numNodePos+(g->quadConn[j][0]+1),
           (i+1)*g->numNodePos+(g->quadConn[j][1]+1),
           (i+1)*g->numNodePos+(g->quadConn[j][2]+1),
           (i+1)*g->numNodePos+(g->quadConn[j][3]+1));
         }

      }
      fclose(fptr);
   }


}
// ##################################################################
// END OF FILE
// ##################################################################
