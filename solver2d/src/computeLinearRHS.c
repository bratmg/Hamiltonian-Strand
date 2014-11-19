// ##################################################################
//
// computeLinearRHS.c
//
// ##################################################################
#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void computeLinearRHS(GRID *g,SOLN *s,double cflnum,double *l2rho)
{
  //
  int    i,j,k,m,f,n;
  int    f1,f2;
  int    is,ie;
  int    iface;
  int    idv;
  int    chainSize;
  int    node1,node2,leftCell,rightCell,icell;
  double ds[2];
  double leftState[NVAR];
  double rightState[NVAR];
  double dleftState[NVAR];
  double drightState[NVAR];
  double consVar[NVAR];
  double dqvar[NVAR];
  double flux[NVAR];
  double gm1     = gamm-1.0;
  double gamma1  = gamm;
  double specRadius;
  double faceVel = 0.;
  double dsnorm,nynx,nx2ny;
  double dtfac;
  double rhoi;
  double x1,y1,x2,y2;  
  double pp;
  double th,qt,eps;
  double dscheck[2];

  int nghost,order;

  // set number of ghost cells
  order = g->order;

                 nghost = 2;
  if(order == 5) nghost = 3;
  //
  // add diagonal term to the linear residual
  //
  for(i=0;i<NVAR*g->ncells;i++) s->r[i]=(s->r0[i]-s->dq[i]);
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  for(i=0;i<g->nchains;i++)
  {
    f1 = g->faceStartPerChain[i];
    f2 = g->faceStartPerChain[i+1];
    m  = nghost;

    for(f=f1;f<f2;f++)
    {
      iface       = g->chainConn[f];
      g->cindx[m] = g->faces[6*iface+2];
      m++;
    }
    //
    // add buffer cells to the chain
    //
    if (g->chainConn[f1]==g->chainConn[f2-1])
    {
      //
      // this is a closed chain
      // make it periodic
      //
      f           = f1+1;
      iface       = g->chainConn[f];
      g->cindx[m] = g->faces[6*iface+2];
      m++;
      chainSize = m;
      
      m = 0;
      for(f=f2-nghost-1;f<f2-1;f++)
      {
        iface=g->chainConn[f];
        g->cindx[m]=g->faces[6*iface+2];
        m++;
      }
    }
    else
    {
      //
      // this is a open chain
      // -ve index indicates necessity to create
      // ghost cells
      //
      if(g->test == 0)
      { 
        if(order==5) 
        {
          m--;
          g->cindx[m]=-g->cindx[m];
          m++;
          g->cindx[m]=-g->cindx[m-3];
          m++;
          g->cindx[m]=-g->cindx[m-5];

          chainSize=m+1;
          m=0;
          g->cindx[m]=-g->cindx[m+5];
          m=1;
          g->cindx[m]=-g->cindx[m+3];
          m=2;
          g->cindx[m]=-g->cindx[m+1];
        }
        else
        {
          m--;
          g->cindx[m]=-g->cindx[m];
          m++;
          g->cindx[m]=-g->cindx[m-3];
          chainSize=m+1;
          m=0;
          g->cindx[m]=-g->cindx[m+3];
          m=1;
          g->cindx[m]=-g->cindx[m+1];
        }
      } //test =0 

      if(g->test==1)
      {
        apply_periodic(&g[0],f1,f2,m);
        chainSize = m+1;
        if(order==5) chainSize = m+2;
      }
    }//open loops

    for(j=0;j<chainSize;j++)
    {
      icell = g->cindx[j];
      if (icell >=0) 
      {
        m = NVAR*icell;
        for(k=0;k<NVAR;k++)
        {
          consVar[k] = s->q[m];
          dqvar[k]   = s->dq[m];
          m++;
        }
        rhoi=1./consVar[0];
        //
        // collect primitive variables in 
        //  
        g->f[j][0]=consVar[0];
        g->f[j][1]=consVar[1];
        g->f[j][2]=consVar[2];
        g->f[j][3]=consVar[3];
        // 
        g->df[j][0]=dqvar[0];
        g->df[j][1]=dqvar[1];
        g->df[j][2]=dqvar[2];
        g->df[j][3]=dqvar[3];
      }
      else
      {
        //
        // ghost cells
        // based on whether they are on the solid boundary on that
        //
        if (j < nghost) 
          iface=g->chainConn[f1];
        else
          iface=g->chainConn[f2-1];
        
        rightCell=g->faces[6*iface+4];
         
        if (rightCell==-2)  /* this is a face on solid wall */
        {
          node1 = g->faces[6*iface];
          node2 = g->faces[6*iface+1];
          
          x1    = g->x[2*node1];
          y1    = g->x[2*node1+1];
          x2    = g->x[2*node2];
          y2    = g->x[2*node2+1];
          
          ds[0] =  (y2-y1);
          ds[1] = -(x2-x1);

          icell =-icell;

          m = NVAR*icell;
          for(k=0;k<NVAR;k++)
          {
            consVar[k] = s->q[m];
            dqvar[k]   = s->dq[m];
            m++;
          }

          dsnorm = ds[0]*ds[0]+ds[1]*ds[1];
          nynx   = ds[0]*ds[1]/dsnorm;
          nx2ny  = (ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
          rhoi   = 1./consVar[0];
          /*  
          g->f[j][0]=consVar[0];
          g->f[j][1]=(-consVar[1]*nx2ny-2*consVar[2]*nynx)*rhoi;
          g->f[j][2]=(consVar[2]*nx2ny-2*consVar[1]*nynx)*rhoi;
          g->f[j][3]=gm1*(consVar[3]-0.5*(consVar[1]*consVar[1]+consVar[2]*consVar[2])*rhoi);
          */
          g->f[j][0]  = consVar[0];
          g->f[j][1]  = (-consVar[1]*nx2ny-2*consVar[2]*nynx);
          g->f[j][2]  = (consVar[2]*nx2ny-2*consVar[1]*nynx);
          g->f[j][3]  = consVar[3];              

          g->df[j][0] = dqvar[0];
          g->df[j][1] = (-dqvar[1]*nx2ny-2*dqvar[2]*nynx);
          g->df[j][2] = (dqvar[2]*nx2ny-2*dqvar[1]*nynx);
          g->df[j][3] = dqvar[3];            
        }
        else 
        {
          g->f[j][0]  = rinf;
          g->f[j][1]  = rinf*s->uinf;
          g->f[j][2]  = rinf*s->vinf;
          g->f[j][3]  = pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);

          g->df[j][0] = 0.;
          g->df[j][1] = 0.;
          g->df[j][2] = 0.;
          g->df[j][3] = 0.;
        } // if rightcell == -2
      } // if icell >= 0
    } // j loop (chainsize)
     
 
    is  = nghost-1;
    ie  = chainSize-1;
    th  = 1./3;
    qt  = 0.25;
    if (g->order==1) qt=0.0;
    eps = 1e-10;

    muscld_deriv(g->f,g->dql,g->dqr,g->f2,g->df,is,ie,th,qt,eps,chainSize,NVAR);
    
    // if(order==5) weno_deriv(g->f,g->dql,g->dqr,g->df,is,ie,eps,chainSize,NVAR);

    n   = is;
    idv = (g->chainConn[f1]==g->chainConn[f2-1]);

    for(f=f1;f<f2-idv;f++)
    {
      iface     = g->chainConn[f];
      leftCell  = g->faces[6*iface+2];
      rightCell = g->faces[6*iface+4];

      for(m=0;m<NVAR;m++)
      {
        if (f==f2-idv-1 && idv==0) 
        {
          dleftState[m]  = g->dql[n][m];
          drightState[m] = g->dqr[n+1][m];
        }
        else
        {
          dleftState[m]  = g->dqr[n+1][m];
          drightState[m] = g->dql[n][m];
        }
      } // m loop
     
      for(j=0;j<NVAR;j++)
      {
        flux[j]=0;
        
        for(k=0;k<NVAR;k++)
          flux[j]+=(((g->ff[iface]).lmat[j][k]*dleftState[k])+
                    ((g->ff[iface]).rmat[j][k]*drightState[k]));
        
      }
      //
      m     = NVAR*leftCell;
      dtfac = cflnum/s->sigma[leftCell];
      for(j=0;j<NVAR;j++)
      {
        s->r[m]-=(flux[j]*dtfac);
        m++;
      }
      
      // if NOT a boundary cell
      if (rightCell > -1) 
      {
        m     = NVAR*rightCell;
        dtfac = cflnum/s->sigma[rightCell];
        for(j=0;j<NVAR;j++)
        {
          s->r[m]+=(flux[j]*dtfac);
          m++;
        }
      } // if rightCell > -1 
      n++;
    } // f loop (f1 to f2)
  } // i loop (nchains)

  *l2rho = 0.;
  for (i=0;i<g->ncells;i++)
  {
    if ((*l2rho) < fabs(s->r[4*i])) 
    {
      icell=i;
      *l2rho=fabs(s->r[4*i]);
    }
  }

}
// ##################################################################
// END OF FILE
// ##################################################################
