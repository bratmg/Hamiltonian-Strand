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
#include "flushToSurface.h"

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
   
   int    i,j,k,ii,jj,kk;
   int    nLayer,nConstLayer;
   char   c,buffer[100];
   char   nodeBase[250],nodeData[250];
   char   connBase[250],connData[250];
   char   boundBase[250],boundData[250];
   char   normalBase[250],normalData[250];   
   char   hybridBase[250],hybridData[250];
   char   *appendStr, intStr[10];
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
   fscanf(fptr,"if hybrid     = %d\n",&iHybrid);
   fscanf(fptr,"nhybrid grid  = %d\n",&nHGrid);
   fscanf(fptr,"iHybPeriodic  = %d\n",&iHybPeriodic);
   fscanf(fptr,"node file     = %s\n",nodeBase); 
   fscanf(fptr,"conn file     = %s\n",connBase);
   fscanf(fptr,"normal file   = %s\n",normalBase);
   fscanf(fptr,"boundary file = %s\n",boundBase);
   fscanf(fptr,"hybrid file   = %s\n",hybridBase);
   fscanf(fptr,"flush surface = %s\n",surfaceType);
   fscanf(fptr,"quad level    = %d\n",&qLevel);
   fscanf(fptr,"num smooth    = %d\n",&numSmooth);
   fscanf(fptr,"smooth type   = %s\n",smoothTechnique);
   fscanf(fptr,"if strands    = %d\n",&iStrand);
   fscanf(fptr,"strand layers = %d\n",&nLayer);
   fscanf(fptr,"const layers  = %d\n",&nConstLayer);
   fscanf(fptr,"init spacing  = %lf\n",&initLen);
   fscanf(fptr,"mesh growth   = %lf\n",&growth);
   fclose(fptr);

   if(iStrand && nLayer < 2) nLayer = 2;

   // load robin fuselage data if fuselage type == robin
   if(strcmp(surfaceType,"robin")==0) 
   {  
      createRobinData();
   }


   // initialize the number of grids
   g = (GRID *) malloc(sizeof(GRID)*nGrid);

   appendStr = ".dat";

   if(nGrid >= 1000)
   {
      printf("#meshgen: Too many subdomains. Need to recode \n");
      printf("#meshgen: the numbering scheme. Stopping. \n");
      exit(1);
   }
   // ===============================================================
   // Loop through number of grids
   // ===============================================================
   for (i = 0; i < nGrid; i++) // loop through total number of grids
   {

      strcpy(nodeData,nodeBase);
      strcpy(connData,connBase);
      strcpy(boundData,boundBase);
      strcpy(normalData,normalBase);
      // Obtain the coord and conn file name if more than one domain
      // i.e., if multiple domains or sub-domains
      // if(1)
      if(nGrid > 1) // need to append 001.dat, 002.dat, etc.
      {
         if(i < 10)
            sprintf(intStr,"00%d",i);
         else if(i < 100)
            sprintf(intStr,"0%d",i);
         else
            sprintf(intStr,"%d",i);

         strcat(nodeData,intStr);
         strcat(nodeData,appendStr);
         strcat(connData,intStr);
         strcat(connData,appendStr);
         strcat(boundData,intStr);
         strcat(boundData,appendStr);
         strcat(normalData,intStr);
         strcat(normalData,appendStr);
      }
      else // single grid (need to append the .dat)
      {
         strcat(nodeData,appendStr);
         strcat(connData,appendStr);
         strcat(normalData,appendStr);
      }

      // ============================================================
      // Read coordinate data
      // ============================================================
      if( (fptr = fopen(nodeData,"r")) == NULL)
      {
         printf("Error: Node data file '%s' missing. \n",nodeData);
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
         printf("Error: Connectivity data file '%s' missing. \n",connData);
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
      // Read boundary data
      // ============================================================
      if( nGrid > 1)
      { 
         if ((fptr = fopen(boundData,"r")) == NULL)
         {
            printf("Error: Boundary data file '%s' missing. \n",boundData);
            exit(1);
         }

         // allocate
         fscanf(fptr,"%d\n",&g[i].boundaryCount);
         g[i].boundaryID = (int *) malloc(sizeof(int)*g[i].boundaryCount);

         for (j = 0; j < g[i].boundaryCount; j++)
            fscanf(fptr,"%d\n",&g[i].boundaryID[j]);

         fclose(fptr);
      }

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


      // ============================================================ 
      // Read node normal data
      // ============================================================
      // allocate
      double dummy;
      int   halfVertm1      = 0.5*(g[i].nVertPerSide-1);
      
      k = g[i].numTriNode + g[i].numTriangle*(3*g[i].nVertPerSide +
                                          3*halfVertm1 + 1 + 
                                          3*halfVertm1*halfVertm1);
      k*=3;

      // allocation an upper bound!
      g[i].nodeNormal = (double *) malloc(sizeof(double)*k);
      for (j = 0; j < k; j++)
         g[i].nodeNormal[j] = -100.;

      if (strcmp(surfaceType,"robin")==0)
      {
         if ((fptr = fopen(normalData,"r")) == NULL)
         {
            printf("Error: Normal data file '%s' missing. \n",normalData);
            exit(1);
         }

         fscanf(fptr,"%d\n",&dummy);

         k = 0;
         for (j = 0; j < 3*g[i].numTriNode; j++)
         {
            fscanf(fptr,"%lf %lf %lf\n",&g[i].nodeNormal[k]
               ,&g[i].nodeNormal[k+1],&g[i].nodeNormal[k+2]);
            k+=3;
         }

         fclose(fptr);
      }


   } // end i for number of grid

   // ============================================================ 
   // Read the structured body fitted data if hybrid exits
   // ============================================================
   if(iHybrid)
   {
      //
      // variable declarations local to this section of the code
      //
      int tempint,indx;

      // error check for more than 1 hybrid grid
      if(nHGrid > 1)
      {
         printf("Error: More than 1 hybrid grid chosen.\n");
         exit(1);
      }

      // initialize hybrid grid structure
      h = (HGRID *) malloc(sizeof(HGRID)*nHGrid);

      //
      // Loop through the total number of hybrid grids
      // and extract position information from data files
      //
      for (i = 0; i < nHGrid; i++)
      {

         // create the filename
         strcpy(hybridData,hybridBase);
         strcat(hybridData,appendStr);

         if( (fptr = fopen(hybridData,"r"))==NULL)
         {
            printf("Error: Hybrid data file '%s' missing. \n",hybridData);
            exit(1);
         }

         // read the data from the files
         fscanf(fptr,"%d\n",&h[i].nPsi);
         fscanf(fptr,"%d\n",&h[i].nEta);
         
         // allocate nodePos
         tempint      = (4*qLevel*(h[i].nPsi-1) + 1)*h[i].nEta;
         h[i].nodePos = (double *) malloc(sizeof(double)*3*tempint);
         for (ii = 0; ii < 3*tempint; ii++)
            h[i].nodePos[ii] = -0.;
         
         // load data into nodepos         
         for (ii = 0; ii < h[i].nEta; ii++)
         {  
            for (jj = 0; jj < h[i].nPsi; jj++)
            {
               indx = (4*qLevel*(h[i].nPsi-1)+1)*ii + jj*(4*qLevel);

               fscanf(fptr,"%lf",&h[i].nodePos[3*indx  ]);
               fscanf(fptr,"%lf",&h[i].nodePos[3*indx+1]);


            } // jj loop
         } // ii loop

         fclose(fptr);

      } // i loop (nHGrid)


   } // if hybrid

}

// ##################################################################
//
// writeTecplot
// reads from input.meshgen and from the node and connectivity data
//
// ##################################################################
void writeTecplot(int gridID, GRID *g)
{

   int      i,j,k,kk,temp;
   int      ii,i1,i2,jj;
   int      ie,iA,iB,iC,iD;
   int      count1,count2,lenCount;
   int     *edgeConn = (int *) malloc(sizeof(int)*2*g->numQuadEdge);
   char     filename[30],fileID[250],intStr[10]; 
   double   dist;
   double  *edgePts  = (double *)  malloc(sizeof(double)*3*g->quadLoop->totLen);   
   double **lenLoop  = (double **) malloc(sizeof(double *)*pow2*g->numTriNode);
   FILE    *fptr;


   // ===============================================================
   // TRIANGLES.DAT
   // ===============================================================   
   strcpy(fileID,folderName);
   strcat(fileID,"triangles.tecdat");
   fptr = fopen(fileID,"w");
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
   // normalize the normals
   // double sum;
   // sum = 0.;
   // k = 0;
   // for (j = 0; j < g->numNodePos; j++)
   // {
   //    sum += g->nodeNormal[3*j  ]*g->nodeNormal[3*j  ] 
   //         + g->nodeNormal[3*j+1]*g->nodeNormal[3*j+1] 
   //         + g->nodeNormal[3*j+2]*g->nodeNormal[3*j+2];
   //    sum = 1./sqrt(sum);

   //    g->nodeNormal[3*j  ]*=sum;
   //    g->nodeNormal[3*j+1]*=sum;
   //    g->nodeNormal[3*j+2]*=sum;

   // }



   strcpy(fileID,folderName);
   strcat(fileID,"quads.tecdat");
   fptr = fopen(fileID,"w");
   fprintf(fptr,"VARIABLES=\"X\",\"Y\",\"Z\" \n");
   fprintf(fptr,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",
      g->numNodePos,g->numQuadConn);

   // write the position information
   k = 0;
   for (j = 0; j <  g->numNodePos; j++)
   {
      fprintf(fptr, "%lf %lf %lf \n", 
         g->allNodePos[k],g->allNodePos[k+1],g->allNodePos[k+2]);
      // fprintf(fptr, "%lf %lf %lf \n", 
      //   g->allNodePos[k  ]+0.1*g->nodeNormal[k],
      //   g->allNodePos[k+1]+0.1*g->nodeNormal[k+1],
      //   g->allNodePos[k+2]+0.1*g->nodeNormal[k+2]);
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
   strcpy(fileID,folderName);
   strcat(fileID,"coord.dat");
   fptr = fopen(fileID,"w");
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
      strcpy(fileID,folderName);
      strcat(fileID,"conn.dat");
      fptr = fopen(fileID,"w");
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
      strcpy(fileID,folderName);
      strcat(fileID,"conn.dat");
      fptr = fopen(fileID,"w");
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
      strcpy(fileID,folderName);
      strcat(fileID,"ofaces.dat");
      fptr = fopen(fileID,"w");
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
      strcpy(fileID,folderName);
      strcat(fileID,"qedges.dat");
      fptr = fopen(fileID,"w");
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
   strcpy(fileID,folderName);
   strcat(fileID,"qloops.dat");
   fptr = fopen(fileID,"w");
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
   strcpy(fileID,folderName);
   strcat(fileID,"iqloops.dat");
   fptr = fopen(fileID,"w");
   if(iStrand) 
   {
      // understand the reason for +2
      temp = (pow2*g->numTriNode+1)*(g->numStrandLayer-1)+g->numQuadConn;

      // include body fitted loops also
      if(iHybrid) temp += (g->nEta-1)*(g->numStrandLayer-1);


      fprintf(fptr,"%d\n",temp);
      for (j = 0; j < temp; j++)
         fprintf(fptr, "%d\n", g->cellLoop->index[j]);
   }
   else
   {
      // temp = pow2*g->numTriNode+1;
      temp = g->niqloops;
      fprintf(fptr,"%d\n",temp);
         
      for (j = 0; j < temp; j++)
         fprintf(fptr, "%d\n", g->quadLoop->index[j]);
         
   }

   
   fclose(fptr);

   // ===============================================================
   // NCOLORS.DAT
   // ===============================================================
   strcpy(fileID,folderName);
   strcat(fileID,"ncolors.dat");
   fptr = fopen(fileID,"w");
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
      int ktemp,stride,offset;
      temp  = pow2*g->numQuadEdge*(g->numStrandLayer-1)
                + g->numStrandLayer*g->numQuadConn + 1;
      int    *faceConn = (int *)     malloc(sizeof(int)*temp);
      double *facePts  = (double *)  malloc(sizeof(double)*3*g->cellLoop->totLen);   

      k      = 0;
      ktemp  = 0;

      stride = pow2*g->numTriNode;
      offset = 0;
      // if(iHybrid) stride += (g->nEta-1);
      //
      // loop through all the colours
      //
      for (i = 0; i <= g->maxCol; i++)      
      {

         sprintf(filename,"loops%d.tecdat", i+1); 
         strcpy(fileID,folderName);
         strcat(fileID,filename);
         fptr = fopen(fileID,"w");

         count1 = 0;
         count2 = 0;
         kk     = 0;
         temp   = pow2*g->numCol[i];  
         if(iHybrid)
         {
            if(i == g->maxCol)
            {
               temp   = g->numCol[i];
               stride = g->nEta-1;
               offset = (g->numStrandLayer-1)*pow2*g->numTriNode;
               ktemp  = 0;
            }
         }
         trace(temp);
         //
         // loop through all the strand layers
         //
         for (jj = 0; jj < g->numStrandLayer-1; jj++)         
         {
            
            k = ktemp;            
            for (j = 0; j < temp; j++)
            {         
               i1 = g->cellLoop->index[offset + k   + jj*(stride+1)];
               i2 = g->cellLoop->index[offset + k+1 + jj*(stride+1)];                                 
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
      strcpy(fileID,folderName);
      strcat(fileID,"nstrands.dat");
      fptr = fopen(fileID,"w");
      fprintf(fptr,"%d\n",g->numStrandLayer);
      fclose(fptr);

      // STRANDS.DAT
      printf("#meshgen: Preparing Tecplot file for strands ...\n");         
      strcpy(fileID,folderName);
      strcat(fileID,"strands.tecdat");
      fptr = fopen(fileID,"w");
      
      k      = (pow2*g->numTriNode+1)*(g->numStrandLayer-1)-1;
      if(iHybrid) k += (g->nEta-1)*(g->numStrandLayer-1);

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
         sprintf(filename,"loops%d.tecdat", i+1);      

         strcpy(fileID,folderName);
         strcat(fileID,filename);
         fptr = fopen(fileID,"w");

         count1 = 0;
         count2 = 0;
         kk     = 0;
         temp   = pow2*g->numCol[i];  
         
         if(iHybrid)
         {
            if(i == g->maxCol)
            {
               temp = g->numCol[i];
            }
         }

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
      strcpy(fileID,folderName);
      strcat(fileID,"brick.tecdat");
      fptr = fopen(fileID,"w");

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



   // ===============================================================
   // domainEdges
   // ===============================================================
   if (nGrid > 1)
   {
      strcpy(fileID,folderName);
      strcat(fileID,"domainEdges.tecdat");
      fptr = fopen(fileID,"w");
      fprintf(fptr,"VARIABLES=\"X\",\"Y\",\"Z\" \n");
      fprintf(fptr,"ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG\n",
      g->numTriNode,g->boundaryCount);

      for (j = 0; j < g->numTriNode; j++)
      {
         fprintf(fptr,"%lf %lf %lf \n",
            g->nodePosTri[3*j],g->nodePosTri[3*j+1],g->nodePosTri[3*j+2]);
      }

      for (j = 0; j < g->boundaryCount; j++)
      {
         fprintf(fptr,"%d %d \n",
            g->triEdge[g->boundaryID[j]][0]+1,g->triEdge[g->boundaryID[j]][1]+1);
      }

      fclose(fptr);
   }

}

// ==================================================================
//
// createOutputFolder.c
//
// creates the name of the subdomain folder and the variable is
// stored in globalVariables.c
// ==================================================================
void createOutputFolder(int gridID)
{
	char intStr[10],command[250];

   // output to screen
   printf("\n====================================================================\n");
   printf("\nDOMAIN %d\n\n",gridID );
   printf("====================================================================\n");
   
   // create folder
   if (gridID < 10)
      sprintf(intStr,"00%d/",gridID);                                                                             
   else if (gridID < 100)
      sprintf(intStr,"0%d/",gridID);
   else if (gridID < 1000)
      sprintf(intStr,"%d/",gridID);

	strcpy(folderName,"./output/domain");
	strcat(folderName,intStr);
	strcpy(command,"mkdir ");
	strcat(command,folderName);
	system(command);

}


// ##################################################################
// Write array to file (utility functions)
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
// END OF FILE
// ##################################################################
