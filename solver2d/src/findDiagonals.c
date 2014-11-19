// ##################################################################
//
// findDiagonals.c
//
// utility function for DADI.c
// ##################################################################
#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define NQ 4

void findDiagonals(GRID *g, SOLN *s,double cflnum,double dt)
{
  int i,j,k,l,m,f,n;
  int mm1,f1,f2,iface,idv,chainSize;
  double ds[2];
  double lmat[NQ][NQ];
  double rmat[NQ][NQ];
  double dsnorm,nynx,nx2ny;
  int node1,node2,leftCell,rightCell,icell;
  double x1,y1,x2,y2,r01,r02,a1,a2,b1,b2,pp,dtfac;
  double linearl2rho;
  double damping;
  //
  // set diagonal term to zero
  //
  for(i=0;i<g->ncells;i++)
  { 
    s->itag[i]=0;
    for(j=0;j<NQ;j++)
      for(k=0;k<NQ;k++)
   
    s->D[i][j][k]=0;
  }
  //
  // first collect the diagonal terms
  //
  for(i=0;i<g->nchains;i++)
  {
    //
    // collect cells on the loop
    // 
    f1  = g->faceStartPerChain[i];
    f2  = g->faceStartPerChain[i+1];
    idv = (g->chainConn[f1]==g->chainConn[f2-1]);
    
    m         = 0;
    chainSize = (f2-idv-f1);
    for(f=f1;f<f2-idv;f++)
    {
      iface     = g->chainConn[f];
      leftCell  = g->faces[6*iface+2];
      rightCell = g->faces[6*iface+4];
      //
      // transpose since jac_roe_ was in f90 ordering
      // remove when jac_roe goes to C
      //
      for(j=0;j<NQ;j++)
        for(k=0;k<NQ;k++)
        {
          lmat[j][k] = (g->ff[iface]).lmat[j][k];
          rmat[j][k] = (g->ff[iface]).rmat[j][k];
        }
      
      for(j=0;j<NQ;j++)
        for(k=0;k<NQ;k++)
          s->D[leftCell][j][k]+=(lmat[j][k]);
      
      // if NOT boundary cell
      if (rightCell > -1)
      {
        for(j=0;j<NQ;j++)
          for(k=0;k<NQ;k++)
            s->D[rightCell][j][k]-=(rmat[j][k]);
      }
    } // for loop (f1 to f2)
  } // i loop (nchains)

  damping = 1.0 + 0.1*cflnum*(g->msweep > 1);
  for(i=0;i<g->ncells;i++)
  {
    dtfac=cflnum/s->sigma[i];

    for(j=0;j<NQ;j++)
      for(k=0;k<NQ;k++)
        s->D[i][j][k]*=dtfac;

    for(j=0;j<NQ;j++)
      s->D[i][j][j]+=damping;
  } // i loop
  
}
// ##################################################################
// END OF FILE
// ##################################################################
