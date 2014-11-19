// ##################################################################
//
// DADI.c
//
// ##################################################################
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define NQ 4

void DADI(GRID *g,SOLN *s,double cflnum,double dt)
{
  //
  int    i,j,k,l,m,f,n;
  int    mm1,f1,f2,iface,idv,chainSize;
  int    isweep;
  int    node1,node2,leftCell,rightCell,icell;
  double ds[2];
  double lmat[NQ][NQ];
  double rmat[NQ][NQ];
  double dsnorm,nynx,nx2ny;
  double x1,y1,x2,y2,r01,r02,a1,a2,b1,b2,pp,dtfac;
  double linearl2rho;
  double damping;
  double tval[NQ];
  //
  // store residuals in case of linear sweep
  //
  for(i=0;i<g->ncells;i++)
  {
    dtfac=cflnum/s->sigma[i];
    for(m=0;m<NVAR;m++)
    {
      s->r[NVAR*i+m]*=dtfac;
      s->r0[NVAR*i+m] = s->r[NVAR*i+m];
      s->dq[NVAR*i+m] = 0.;
    }
  } // i loop
  //
  // find diagonal matrix terms
  //
  findDiagonals(g,s,cflnum,dt);
  for (isweep=0;isweep<g->msweep;isweep++)
  {
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
        // collect left and right matrices for state
        //
        for(j=0;j<NQ;j++)
          for(k=0;k<NQ;k++)
          {
            lmat[j][k] = (g->ff[iface]).lmat[j][k];
            rmat[j][k] = (g->ff[iface]).rmat[j][k];
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
              g->A[m][j][k]+=(rmat[j][k]);
              g->C[mm1][j][k]-=(lmat[j][k]);
            }
            g->F[m][j]=s->r[NVAR*leftCell+j];
          } // j loop
        }
        //
        // open loop
        //
        else
        {
          mm1=(m-1);
          if (rightCell < 0 && idv==0 && f==f2-1) m--;
        
          for(j=0;j<NQ;j++)
          {
            for(k=0;k<NQ;k++)
            {
              if (mm1 > -1 && rightCell > -1)
              {
                g->A[m][j][k]+=(rmat[j][k]);
                g->C[mm1][j][k]-=(lmat[j][k]);
              }
            } // k loop
            g->F[m][j]=s->r[NVAR*leftCell+j];
          } // j loop
        } // if (idv==1)  
        
        m++;
      } // f loop (f1 to f2)

      m         = 0;
      chainSize = chainSize-(idv==0);
      for(f=f1;f<f2-1;f++)
      {
        iface    = g->chainConn[f];
        leftCell = g->faces[6*iface+2];
        dtfac    = cflnum/s->sigma[leftCell];
        for(j=0;j<NQ;j++)
          for(k=0;k<NQ;k++)
          {
            g->A[m][j][k]*=dtfac;
            g->B[m][j][k]=s->D[leftCell][j][k];
            g->C[m][j][k]*=dtfac;
          }
          
        m++;
      } //  f loop (f1 to f2)
      //
      // invert using appropriate banded block solver
      //
      if (idv==1) blockTridagPeriodic4(g->A,g->B,g->C,g->F,chainSize,NQ);
      if (idv==0) blockTridag4(g->A,g->B,g->C,g->F,chainSize,NQ);
      //
      // multiply by D block
      // reassign values back at the unknown locations
      //
      m=0;
      for(f=f1;f<f2-1;f++)
      {
        iface    = g->chainConn[f];
        leftCell = g->faces[6*iface+2];
        if (s->itag[leftCell] < 1) 
        {
          for(j=0;j<NQ;j++) tval[j]=0;
          
          for(j=0;j<NQ;j++)
            for(k=0;k<NQ;k++)
              tval[j]+=s->D[leftCell][j][k]*g->F[m][k];
        
          s->itag[leftCell]++;
          for(j=0;j<NQ;j++)
            s->r[NVAR*leftCell+j]=tval[j];
        }
        else
        {
          for(j=0;j<NQ;j++)
            s->r[NVAR*leftCell+j]=g->F[m][j];
        }
        m++;
      } // f loop (f1 to f2)            
    } // i sweep

    for(i=0;i<g->ncells*NVAR;i++) s->dq[i]+=s->r[i];
    
    computeLinearRHS(g,s,cflnum,&linearl2rho);
    if (g->msweep > 1) tracef(linearl2rho);
  }
  //
  // update q
  //
  m=0;
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
