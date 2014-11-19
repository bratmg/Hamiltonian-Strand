// ##################################################################
//
// preprocess.c
//
// Preprocessing code
// Written by Dr. Jayanarayanan Sitaraman
// Note: When accessing elements of 'x', the access is random.
//       Efforts can be focussed towards a more linear accessing
//       of data to better utilize cache and increase efficiency
//       of parallelization
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define nearBodyRadius 3.0 // hardcoded as 3.0


// ==================================================================  
// Preprcoess the inputs obtained in reagrid 
// (i.e., matlab mesh output)
// ==================================================================
void preprocess(GRID *g)
{
  int    i,j,k,m,i1,i2;
  int    leftCell,rightCell;
  int    leftFaceIndx,rightFaceIndx;
  int    f,f1,f2,ff1,ff2;
  int    icell,iface;
  double r1,r2;
  double x1,x2,y1,y2,dvedge;
  FILE   *fp;
  double volmin,volmax;
  int tmp,cell1,cell2,cell3,cell4;
  double tmpn;
  double bf1n[5000],bf2n[5000],bf3n[5000],bf4n[5000];
  //double bf1c[500],bf2c[500],bf3c[500],bf4c[500];



  // Perform smoothing of the initial grid
  //if (!g->visc) smoothGrid(g,200);

// ==================================================================
// Begin commented section by Jayanarayanan Sitaraman
// ==================================================================  
  //find_faces(g->conn,&(g->neig),&(g->faces),&(g->nfaces),g->ncells,4);
  /*
  fp=fopen("QuadData/qedges1.dat","w");
  fprintf(fp,"%d\n",(g->nfaces));
  for(i=0;i<g->nfaces;i++)
    fprintf(fp,"%d %d %d %d %d %d\n",(g->faces[6*i]),
     (g->faces[6*i+1]),
     (g->faces[6*i+2]),
     (g->faces[6*i+4]),
     (g->faces[6*i+3]),
     (g->faces[6*i+5]));
  fclose(fp);
  */
  /*
  for(i=0;i<g->nfaces;i++)
    {
      leftCell=g->faces[6*i+2];
      leftFaceIndx=g->faces[6*i+3];
  */  
  //trace(g->nfaces);
  //
  // find boundary faces
  // this code-snippet is not general, its 
  // just for external flow problems where the body is constrained
  // within a radius of nearBodyRadius from origin
// ==================================================================
// End commented section by Jayanarayanan Sitaraman
// ==================================================================  

   
  if(g->test==1) periodic_bc(&g[0]);

  if(g->test==0) 
  {
    // Initialize the boundary faces
    g->nbfaces=0;
    // Loop through the total number of faces
    for (i=0;i<g->nfaces;i++)
    {
       // recall swap for left face = -1 was performed in readGrid
       if (g->faces[6*i+4] == -1) 
       {
          i1 = g->faces[6*i];
          i2 = g->faces[6*i+1];
          r1 = sqrt(g->x[2*i1]*g->x[2*i1] + g->x[2*i1+1]*g->x[2*i1+1]);
          r2 = sqrt(g->x[2*i2]*g->x[2*i2] + g->x[2*i2+1]*g->x[2*i2+1]);

          // if average radius of the edge face is less that near
          // body radius, then increase boundary face count
          if (0.5*(r1+r2) < nearBodyRadius) g->nbfaces++;
       }
    }
  } //test=0

  // Initializing pointers (list of boundary faces, cell volume
  // list of neighbours, cell to face, cell to chain )
  
  if(g->test==0) g->bfaces  = (int *)malloc(sizeof(int)*g->nbfaces);
  g->vol     = (double *) malloc(sizeof(double)*g->ncells);
  g->neig    = (int *)    malloc(4*sizeof(int)*g->ncells);
  g->c2f     = (int *)    malloc(4*sizeof(int)*g->ncells);
  g->c2chain = (int *)    malloc(3*sizeof(int)*g->ncells);

  // initialize the cell volume for all the cells (quadilateral)
  for(i=0;i<g->ncells;i++) g->vol[i]=0.;

  // loop through the number of faces (from Qedges)
  m=0;
  for (i=0;i<g->nfaces;i++)
  {
      leftCell  = g->faces[6*i+2];
      rightCell = g->faces[6*i+4];
      leftFaceIndx  = g->faces[6*i+3];
      rightFaceIndx = g->faces[6*i+5];
      
      // Make a list of faces for a given cell
      g->c2f[4*leftCell + leftFaceIndx] = i;

      if (rightCell > -1) 
      {
         g->c2f[4*rightCell+rightFaceIndx] = i;
      } 

      // collect neighbors
      // Make a list of neighbours  for a given cell
      g->neig[4*leftCell+leftFaceIndx] = rightCell;
      
      if (rightCell > -1) 
      {
         g->neig[4*rightCell+rightFaceIndx] = leftCell;
      }

      //
      // indices of edge nodes
      //
      i1 = g->faces[6*i];
      i2 = g->faces[6*i+1];

      // corresponding x and y values
      x1 = g->x[2*i1];
      y1 = g->x[2*i1+1];
      x2 = g->x[2*i2];
      y2 = g->x[2*i2+1];

      dvedge = (x1*y2-x2*y1)*0.5;
      
      //
      // compute cell volumes
      //
      g->vol[leftCell] += dvedge;

      if (rightCell > -1) g->vol[rightCell] -= dvedge;

      if(g->test==0)
      { // if a boundary face
         if (g->faces[6*i+4] == -1) 
         {
           i1 = g->faces[6*i];
           i2 = g->faces[6*i+1];

           r1 = sqrt(g->x[2*i1]*g->x[2*i1] + g->x[2*i1+1]*g->x[2*i1+1]);
           r2 = sqrt(g->x[2*i2]*g->x[2*i2] + g->x[2*i2+1]*g->x[2*i2+1]);

           if (0.5*(r1+r2) < nearBodyRadius) 
           {               
              g->bfaces[m] = i;
              m++; // running counter update

              g->faces[6*i+4]        = -(g->visc+2);
              icell                  = g->faces[6*i+2]; // used?
              iface                  = g->faces[6*i+3]; // used?
              g->neig[3*icell+iface] = -(g->visc+2);
           }
         }
      } //test==0
   }//nface

   // Initial values for volmin and volmax  
   volmin = 1E15;
   volmax = 0.;

   // evaluate the minimum and maximum cell volume
   for (i=0;i<g->ncells;i++)
   {
      volmin = (volmin < g->vol[i]) ? volmin : g->vol[i]; 
      volmax = (volmax > g->vol[i]) ? volmax : g->vol[i];
   }

  tracef(volmin);
  tracef(volmax);
  tracef(g->vol[0]);


// ==================================================================
// find cell to chain connectivity (cell to loops)
// ==================================================================
   // loop through all the cells and initialialize c2chain
   for(i=0;i<g->ncells;i++) 
   { 
      g->c2chain[2*i]   = -1;
      g->c2chain[2*i+1] = -1;
   }

   // loop through all the chains (inner and outer loops)
   for(i=0;i<g->nchains;i++)
   {
      // Starting and ending index of a chain 
      f1 = g->faceStartPerChain[i];
      f2 = g->faceStartPerChain[i+1];

      // loop through the chain indices
      for(f=f1;f<f2-1;f++)
      {
         iface    = g->chainConn[f]; // index of the face/edge
         leftCell = g->faces[6*iface+2]; // corresponding left cell 

         if (g->c2chain[2*leftCell] > -1)
         {
            g->c2chain[2*leftCell+1] = i;
         }
         else
         {
            g->c2chain[2*leftCell]   = i;
         }
      }
   }

// ==================================================================
// write to verify
// ==================================================================
  if(g->test==0)
  {
    fp = fopen("bface.dat","w");
    for(i=0;i<g->nbfaces;i++)
      {
        iface = g->bfaces[i];
        i1    = g->faces[6*iface];

        fprintf(fp,"%f %f\n",g->x[2*i1],g->x[2*i1+1]);

        i1    = g->faces[6*iface+1];
        
        fprintf(fp,"%f %f\n",g->x[2*i1],g->x[2*i1+1]);

      }
    fclose(fp);
  }
}
// ##################################################################
// END OF FILE
// ##################################################################
