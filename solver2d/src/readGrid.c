// ##################################################################
//
// readGrid.c
//
// Function that reads the grid based on the outputs of the Matlab
// mesh generation code
//
// Data files read from:
//  - coord.dat   (data coordinate points)
//  - conn.dat    (node connectivity information for the triangles)
//  - qedges.dat  (writes the QEdge matrix, 6 cols, refer documentation)
//  - ncolors.dat (number of loops of each colour)
//  - iqloops.dat (index of inner and outer loops of each node)
//  - qloops.dat  (Inner and outer loops  for each node)
//
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <stdlib.h>
#define NQ 4

void readGrid(GRID *g)
{
   FILE   *fp;
   char   c;
   int    i,l1,l2,l3,j;
   int    workSize;
   double fdummy;

// ==================================================================
// read coordinates
// ==================================================================
   fp=fopen("QuadData/coord.dat","r"); 
   fscanf(fp,"%d",&(g->nnodes));
  
   // allocate the required space for the coordinate points
   // for the x and y cartesian position
   g->x=(double *) malloc(sizeof(double)*2*g->nnodes);
   for(i=0;i<g->nnodes;i++)
      fscanf(fp,"%lf %lf %lf\n",&(g->x[2*i]),&(g->x[2*i+1]),&(fdummy));
   fclose(fp);

// ==================================================================
// read connectivity data
// set of three node IDs that make up a quadrilateral
// ==================================================================
   fp=fopen("QuadData/conn.dat","r");
   fscanf(fp,"%d",&(g->ncells));
  
   g->conn=(int *) malloc(sizeof(int)*4*g->ncells);
   for(i=0;i<g->ncells;i++)    
      fscanf(fp,"%d %d %d %d\n",&(g->conn[4*i]),&(g->conn[4*i+1]),
     &(g->conn[4*i+2]),&(g->conn[4*i+3]));
    fclose(fp);

// ==================================================================
// read edges (faces)
// This is the Qedges matrix in the Matlab mesh generation
// NOTE: The ordering of the g->faces is different from the 
//       initial Matlab generation
//
// g->faces: 
//  0 - node ID 1 
//  1 - node ID 2
//  2 - left cell ID 
//  3 - left face index
//  4 - right cell ID 
//  5 - right face index
// ==================================================================
   fp=fopen("QuadData/qedges.dat","r");
   fscanf(fp,"%d",&(g->nfaces));
   g->faces=(int *) malloc(sizeof(int)*6*g->nfaces);
   for(i=0;i<g->nfaces;i++)
   {
      fscanf(fp,"%d %d %d %d %d %d\n",
       &(g->faces[6*i]),    // node ID of one edge
       &(g->faces[6*i+1]),  // node ID of other edge
       &(g->faces[6*i+2]),  // left  cell
       &(g->faces[6*i+4]),  // right cell
       &(g->faces[6*i+3]),  // left  face index (0-3)
       &(g->faces[6*i+5])); // right face index (0-3)

      //
      // swap if left cell = -1 
      // Therefore as a result all -1s are not in [6*i+4]
      // I.E., all boundary edges have their right face to the outside
      if (g->faces[6*i+2] == -1) 
     {
         swap(g->faces[6*i],  g->faces[6*i+1]);
         swap(g->faces[6*i+2],g->faces[6*i+4]);
         swap(g->faces[6*i+3],g->faces[6*i+5]);  
      }
   }
   fclose(fp);

// ==================================================================
// read colours (i.e., number of loops per colour)
// ==================================================================
   fp=fopen("QuadData/ncolors.dat","r");
   fscanf(fp,"%d",&(g->ncolors));
   g->chainsPerColor=(int *) malloc(sizeof(int)*g->ncolors);
   for(i=0;i<g->ncolors;i++)
      fscanf(fp,"%d\n",&(g->chainsPerColor[i]));
   fclose(fp);

// ==================================================================
// read chain information (chain = loops)
// ==================================================================
   fp = fopen("QuadData/iqloops.dat","r");
   fscanf(fp,"%d\n",&(g->nchains));
   
   g->nchains--; // not sure why it is reduced by 1
   g->faceStartPerChain=(int *) malloc(sizeof(int)*(g->nchains+1));
   g->nmaxchain=0;

   for(i=0;i<g->nchains+1;i++)
   {
      fscanf(fp,"%d\n",&(g->faceStartPerChain[i]));
      if (i > 0) 
      {
         g->nmaxchain = max(g->nmaxchain,
         g->faceStartPerChain[i]-g->faceStartPerChain[i-1]);
      }
   }
   fclose(fp);
   //
   trace(g->nmaxchain);

   workSize = g->nmaxchain+6;

// ==================================================================
// Allocate memory
// ==================================================================   
   g->ql    = (double **)  malloc(sizeof(double *)*(workSize));
   g->qr    = (double **)  malloc(sizeof(double *)*(workSize));
   g->dql   = (double **)  malloc(sizeof(double *)*(workSize));
   g->dqr   = (double **)  malloc(sizeof(double *)*(workSize));
   g->f     = (double **)  malloc(sizeof(double *)*(workSize));
   g->fv    = (double **)  malloc(sizeof(double *)*(workSize));
   g->df    = (double **)  malloc(sizeof(double *)*(workSize));
   g->f2    = (double **)  malloc(sizeof(double *)*(workSize));
   g->cindx = (int *)      malloc(sizeof(int)*(workSize));
   g->ctype = (int *)      malloc(sizeof(int)*(workSize));
   g->A     = (double ***) malloc(sizeof(double **)*(workSize));
   g->B     = (double ***) malloc(sizeof(double **)*(workSize));
   g->C     = (double ***) malloc(sizeof(double **)*(workSize));  
   g->F     = (double **)  malloc(sizeof(double *)*(workSize));
   g->Q     = (double **)  malloc(sizeof(double *)*(workSize));
   //
   for (i=0;i<workSize;i++)
   {
      g->ql[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->qr[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->dql[i] = (double *)  malloc(sizeof(double)*NVAR);
      g->dqr[i] = (double *)  malloc(sizeof(double)*NVAR);
      g->f[i]   = (double *)  malloc(sizeof(double)*NVAR);
      g->fv[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->df[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->f2[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->A[i]   = (double **) malloc(sizeof(double *)*NQ);
      g->B[i]   = (double **) malloc(sizeof(double *)*NQ);
      g->C[i]   = (double **) malloc(sizeof(double *)*NQ);
      g->F[i]   = (double *)  malloc(sizeof(double)*NQ);
      g->Q[i]   = (double *)  malloc(sizeof(double)*NQ);
      for(j=0;j<NQ;j++)
     {
        g->A[i][j] = (double *) malloc(sizeof(double)*NQ);
        g->B[i][j] = (double *) malloc(sizeof(double)*NQ);
        g->C[i][j] = (double *) malloc(sizeof(double)*NQ);
     } // j loop
   } // i loop

// ==================================================================
// read chain information (chain connectivity)
// nchainFaces = nloops (from Matlab)
// ================================================================
  fp=fopen("QuadData/qloops.dat","r");
  fscanf(fp,"%d",&(g->nchainFaces));
  g->chainConn=(int *)malloc(sizeof(int) * g->nchainFaces);

  for (i=0;i<g->nchainFaces;i++)
     fscanf(fp,"%d",&(g->chainConn[i]));
  fclose(fp);

  printf("#ham2d: Finished reading files\n");

  trace(g->nnodes);
  trace(g->ncells);
  trace(g->nfaces);
  trace(g->ncolors);
  trace(g->nchains);
  trace(g->nchainFaces);
}
// ##################################################################
// END OF FILE
// ##################################################################
