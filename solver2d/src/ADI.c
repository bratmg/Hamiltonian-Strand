// ##################################################################
//
// ADI.c
//
// ##################################################################
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define NQ 4

void ADI(GRID *g,SOLN *s,double cflnum,double dt)
{
   //
   int    i,j,k,l,m,f,n;
   int    isweep;
   int    mm1,f1,f2,iface,idv,chainSize;
   int    node1,node2,leftCell,rightCell,icell;
   double ds[2];
   double lmat[NQ][NQ];
   double rmat[NQ][NQ];
   double dsnorm,nynx,nx2ny;
   double x1,y1,x2,y2,r01,r02,a1,a2,b1,b2,pp,dtfac;
   double linearl2rho;
   //
   // one loop per chain to evaluate fluxes
   // on all the faces in the chain
   //
   for(i=0;i<g->ncells;i++)
   {
      dtfac=cflnum/s->sigma[i];
      for(m=0;m<NVAR;m++)
      {
         s->r[NVAR*i+m] *= dtfac;
         s->r0[NVAR*i+m] = s->r[NVAR*i+m]; // save old residual
         s->dq[NVAR*i+m] = 0.; // reset dq to 0
      }
   } // i loop (ncells)
 
   // number of times the sweep is performed
   for(isweep=0;isweep<g->msweep;isweep++)
   {
      //
      // Find the A,B,C matrix and F vector for each chain
      //
      for(i=0;i<g->nchains;i++)
      {
         //
         // zero out block matrices
         //
         for(m=0;m<g->nmaxchain+5;m++)
            for(k=0;k<NQ;k++)
            {
               g->F[m][k]=0;
               for(j=0;j<NQ;j++)
                  g->A[m][k][j]=g->B[m][k][j]=g->C[m][k][j]=0;
            }
         //
         // collect cells on the loop
         // 
         f1        = g->faceStartPerChain[i];
         f2        = g->faceStartPerChain[i+1];
         idv       = (g->chainConn[f1]==g->chainConn[f2-1]);
         m         = 0;
         chainSize = (f2-idv-f1); // idv checks open/closed loop
         //
         // make matrices A,B,C and vector F for
         // inversion.
         //
         for(f=f1;f<f2-idv;f++)
         {
            iface     = g->chainConn[f];
            leftCell  = g->faces[6*iface+2];
            rightCell = g->faces[6*iface+4];
            // 
            // construct left and right matrices for this face
            //
            for(j=0;j<NQ;j++)
               for(k=0;k<NQ;k++)
               {
                  lmat[j][k] = (g->ff[iface].lmat[j][k]);
                  rmat[j][k] = (g->ff[iface].rmat[j][k]);
               }
     
            // 
            // closed loop
            //
            if (idv==1) 
            {
               mm1=(m==0)?chainSize-1:m-1;
               for(j=0;j<NQ;j++)
               {
                  for(k=0;k<NQ;k++)
                  {
                     g->B[m][j][k]   += (lmat[j][k]);
                     g->B[mm1][j][k] -= (rmat[j][k]);
                     g->A[m][j][k]   += (rmat[j][k]);
                     g->C[mm1][j][k] -= (lmat[j][k]);
                  }
                  g->F[m][j]=s->r[NVAR*leftCell+j];
               }
            }
            //
            // open loop
            //
            else
            {
               mm1 = (m-1);
               if (rightCell==-2)
               {
                  node1  = g->faces[6*iface];
                  node2  = g->faces[6*iface+1];
                  
                  x1     = g->x[2*node1];
                  y1     = g->x[2*node1+1];
                  x2     = g->x[2*node2];
                  y2     = g->x[2*node2+1];
                  
                  ds[0]  = (y2-y1);
                  ds[1]  = -(x2-x1);
                  dsnorm = ds[0]*ds[0]+ds[1]*ds[1];
                  
                  nynx   = ds[0]*ds[1]/dsnorm;
                  nx2ny  = (ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
                  a1     = -nx2ny;
                  b1     = -2*nynx;
                  a2     = -2*nynx;
                  b2     = nx2ny;
                  for (j=0;j<NQ;j++)
                  {
                     r01        = rmat[j][1];
                     r02        = rmat[j][2];
                     rmat[j][1] = a1*r01+b1*r02;
                     rmat[j][2] = a2*r01+b2*r02;
                  }
               } // rightcell==-2
               //
               if (rightCell < 0 && idv==0 && f==f2-1) m--;
               //
               for(j=0;j<NQ;j++)
               {
                  for(k=0;k<NQ;k++)
                  {
                     g->B[m][j][k] += (lmat[j][k]);
                     if (mm1 > -1 && rightCell > -1)
                     {
                        g->A[m][j][k]   += (rmat[j][k]);
                        g->B[mm1][j][k] -= (rmat[j][k]);
                        g->C[mm1][j][k] -= (lmat[j][k]);
                     }
                     else
                     {
                        if (rightCell==-2) g->B[m][j][k]+=(rmat[j][k]);
                     }
                  }
                  g->F[m][j]=s->r[NVAR*leftCell+j];
               } // for loop (j)
            } // if (idv==1)
            m++;
         } // f loop (f1 to f2)
         
         m         = 0;
         chainSize = chainSize-(idv==0);

         // multiply A,B,C matrix by delta_t 
         // and add 'I' to B matrix
         for(f=f1;f<f2-1;f++)
         {
            iface    = g->chainConn[f];    
            leftCell = g->faces[6*iface+2];
            dtfac    = cflnum/s->sigma[leftCell];
            for(j=0;j<NQ;j++)
               for(k=0;k<NQ;k++)
               {
                  g->A[m][j][k] *= dtfac;
                  g->B[m][j][k] *= dtfac;
                  g->C[m][j][k] *= dtfac;
               }
            
            // add the 'I' matrix
            for(j=0;j<NQ;j++)
               g->B[m][j][j]+=1.0;
               
            m++;
         } // f loop (f1 to f2)
         //
         // invert using appropriate banded block solver
         //
         if (idv==1) blockTridagPeriodic4(g->A,g->B,g->C,g->F,chainSize,NQ);
         if (idv==0) blockTridag4(g->A,g->B,g->C,g->F,chainSize,NQ);
         //
         // reassign values back at the unknown locations
         //
         m=0;
         for(f=f1;f<f2-1;f++)
         {
            iface    = g->chainConn[f];
            leftCell = g->faces[6*iface+2];
            
            for(j=0;j<NQ;j++)
               s->r[NVAR*leftCell+j]=g->F[m][j];
            
            m++;
         }      
      } // i loop (nchains)

      for(i=0;i<g->ncells*NVAR;i++) s->dq[i]+=s->r[i];
      
      computeLinearRHS(g,s,cflnum,&linearl2rho);

      //tracef(linearl2rho);
   } // for loop (isweep)
   //
   // update q
   //
   //tracef(s->dq[0]);
   m = 0;
   for(i=0;i<g->ncells;i++)
      for(j=0;j<NVAR;j++)
      {
         s->q[m]+=s->dq[m];
         m++;
      }
      
}
// ##################################################################
// END OF FILE
// ##################################################################
