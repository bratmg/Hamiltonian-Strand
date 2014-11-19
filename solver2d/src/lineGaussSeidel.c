// ##################################################################
//
// lineGaussSeidel.c
//
// ##################################################################
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define NQ 4

void lineGaussSeidel(GRID *g, SOLN *s, double cflnum, double dt)
{
  //
  int i,j,k,l,m,f,n,ii;
  int mm1,f1,f2,iface,idv,chainSize,idiff;
  int iface1,jface;
  int jleftCell,jrightCell;
  double ds[2];
  double lmat[NQ][NQ];
  double rmat[NQ][NQ];
  double dsnorm,nynx,nx2ny;
  int node1,node2,leftCell,rightCell,icell;
  double x1,y1,x2,y2,r01,r02,a1,a2,b1,b2,pp,dtfac;
  double linearl2rho,l2rho,l2rhodq;
  double damping;
  double tval[NQ];
  int isweep;
  FILE *fp1,*fp2,*fp3,*fp4,*fp5;
  //
  /*   fp1=fopen("A.dat","w"); */
  /*   fp2=fopen("B.dat","w"); */
  /*   fp3=fopen("C.dat","w"); */
  /*   fp4=fopen("F.dat","w"); */
  /*   fp5=fopen("X.dat","w"); */
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  //outputr(g,s);
  for(i=0;i<g->ncells;i++)
  {
    dtfac=cflnum/s->sigma[i];
    for(m=0;m<NVAR;m++)
    {
      s->r[NVAR*i+m]  *= dtfac;
      s->r0[NVAR*i+m]  = s->r[NVAR*i+m];
      s->dq[NVAR*i+m]  = 0.;
    } // m loop
  } // i loop
  //
  // find diagonals
  //
  findDiagonals(g,s,cflnum,dt);
  //
  for (isweep=0;isweep<g->msweep;isweep++)
  {
    for (i=0;i<g->ncells*NVAR;i++) s->ddq[i]=0;

    //
    // loop over chains (2 to nchains-1)
    //
    for (ii=0;ii<2*g->nchains-1;ii++)
    {
      //    outputdq(g,s);
      i     = (ii<g->nchains) ? ii: (2*g->nchains-2-ii);
      idiff = (ii < g->nchains) ? 1: -1;
      //
      // zero out block matrices
      //
      for (m=0;m<g->nmaxchain+6;m++)
        for (k=0;k<NQ;k++)
        {
          g->F[m][k]=0;
          for (j=0;j<NQ;j++)
            g->A[m][k][j]=g->B[m][k][j]=g->C[m][k][j]=0;
        }
      //
      // collect cells on the loop
      // 
      f1 = g->faceStartPerChain[i];
      f2 = g->faceStartPerChain[i+1];
      //
      // here idv == 1 means a full chain 
      //
      idv       = (g->chainConn[f1]==g->chainConn[f2-1]);
      m         = 0;
      chainSize = (f2-idv-f1);

      for(f=f1;f<f2-idv;f++)
      {
        iface     = g->chainConn[f];
        leftCell  = g->faces[6*iface+2];
        rightCell = g->faces[6*iface+4];
        //
        for(j=0;j<NQ;j++)
          for(k=0;k<NQ;k++)
          {
            lmat[j][k]=(g->ff[iface]).lmat[j][k];
            rmat[j][k]=(g->ff[iface]).rmat[j][k];
          }
        //
        // closed chain
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
          }
        }
        //
        // open chain
        //
        else
        {
          mm1 = (m-1);
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
            }
            g->F[m][j]=s->r[NVAR*leftCell+j];
          }
        } // end if (idv)  
        m++;
      } // end for  (f1 to f2)

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
        
        iface1 = g->chainConn[f+1];
        icell  = leftCell;
        for (j=0;j<4;j++)
        {
          jface=g->c2f[4*icell+j];
          if (jface != iface  && jface != iface1)
          {
            jleftCell  = g->faces[6*jface+2];
            jrightCell = g->faces[6*jface+4];
            
            if (icell==jleftCell && jrightCell > -1)
            {
              axb1((g->ff[jface]).rmat,&(s->ddq[NVAR*jrightCell]),
                g->F[m],-dtfac,NVAR);
            }
            else if (icell==jrightCell)
            {
              axb1((g->ff[jface]).lmat,&(s->ddq[NVAR*jleftCell]),
                g->F[m],dtfac,NVAR);
            }
          }
        } // j loop
        m++;
      } // end for  (f1 to f2)
      //
      // invert using appropriate banded block solver
      //
      if (idv==1) blockTridagPeriodic4(g->A,g->B,g->C,g->F,chainSize,NQ);
      if (idv==0) blockTridag4(g->A,g->B,g->C,g->F,chainSize,NQ);
      //
      // set ddq can be to be updated by linear RHS if called
      //
      m = 0;
      for(f=f1;f<f2-1;f++)
      {
        //
        iface = g->chainConn[f];
        icell = g->faces[6*iface+2];
        for(j=0;j<NQ;j++) s->ddq[NVAR*icell+j]=g->F[m][j];
        m++;
      }
    } // for loop (ii, nchains)

    for(i=0;i<g->ncells*NVAR;i++) s->dq[i]+=s->ddq[i];
    
    l2rho = 0.;
    icell = 0;
    for(i=0;i<g->ncells;i++)
    {
      if ((l2rho) < fabs(s->ddq[4*i])) 
      {
        icell = i;
        l2rho = fabs(s->ddq[4*i]);
      }
    }

    // Multiple sweeps 
    if (g->msweep > 1) 
    {
      computeLinearRHS(g,s,cflnum,&linearl2rho);
      tracef(linearl2rho);
    }
      //if (g->msweep > 1) tracef(linearl2rho);
  } // isweep loop

  l2rhodq = 0.;
  for(i=0;i<g->ncells;i++)
  {
    if ((l2rhodq) < fabs(s->dq[4*i])) 
    {
      icell   = i;
      l2rhodq = fabs(s->dq[4*i]);
    }
  }
  //
  // update q
  //
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
