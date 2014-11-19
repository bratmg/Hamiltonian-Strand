// ##################################################################
//
// computeRHSk.c
//
// ##################################################################
#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"


void computeRHSk(GRID *g,SOLN *s,double *l2rho)
{
  //
  int i,j,k,m,f,n;
  int f1,f2;
  int is,ie;
  int iface;
  int idv;
  int chainSize;
  int itype;
  int imode=1;
  double ds[2];
  double leftState[NVAR];
  double rightState[NVAR];
  double leftState0[NVAR];
  double rightState0[NVAR];
  double lmat[NVAR][NVAR];
  double rmat[NVAR][NVAR];
  double consVar[NVAR];
  double flux[NVAR];
  double gm1=gamm-1.0;
  double gamma1=gamm;
  double specRadius;
  double faceVel=0.;
  double dsnorm,nynx,nx2ny;
  double rhoi;
  int node1,node2,leftCell,rightCell,icell;
  double x1,y1,x2,y2;  
  double pp;
  double th,qt,eps;
  double dscheck[2];

  int nghost,order,iflag;

  // set number of ghost cells
  order = g->order;
                nghost = 2;
  if(order ==5) nghost = 3;
  //
  // zero out residual and spectral radii
  //
  dscheck[0] = dscheck[1] = 0;

  for(i=0;i<NVAR*g->ncells;i++) s->r[i]     = 0.0;
  for(i=0;i<g->ncells;i++)      s->sigma[i] = 0.0;
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  for(i=0;i<g->nchains;i++)
  {
      iflag = 0;
      f1    = g->faceStartPerChain[i];
      f2    = g->faceStartPerChain[i+1];

      m = nghost;
      for(f=f1;f<f2;f++)
      {
         iface       = g->chainConn[f];
         g->cindx[m] = g->faces[6*iface+2];
         g->ctype[m] = 1;
         m++;
      }
      //
      // add buffer cells to the chain
      //
      if (g->chainConn[f1]==g->chainConn[f2-1])
      {
         iflag = 0;
         //
         // this is a closed chain
         // make it periodic
         //
         f           = f1+1;
         iface       = g->chainConn[f];
         g->cindx[m] = g->faces[6*iface+2];
         g->ctype[m] = 1;
         m++;

         chainSize = m;
         m         = 0;

         for(f=f2-nghost-1;f<f2-1;f++)
         {
            iface       = g->chainConn[f];
            g->cindx[m] = g->faces[6*iface+2];
            g->ctype[m] = 1;
            m++;
         }
      }
      else
      {
         iflag = 1;
         //construct ghost cell using other side cell index
         //solid bc
         if(g->test == 0){
            if(order==5) //WENO 5
            {     
               m--;    
               g->cindx[m] = -g->cindx[m];
               g->ctype[m] = -1;
               m++;
               g->cindx[m] = -g->cindx[m-3];
               g->ctype[m] = -1;
               m++;
               g->cindx[m] = -g->cindx[m-5];
               g->ctype[m] = -1;
               chainSize = m+1;
               m = 0;
               g->cindx[m] = -g->cindx[m+5];
               g->ctype[m] = -1;
               m = 1;
               g->cindx[m] = -g->cindx[m+3];
               g->ctype[m] = -1;
               m = 2;
               g->cindx[m] = -g->cindx[m+1];
               g->ctype[m] = -1;
            }
            else
            {
               m--;
               g->cindx[m] = -g->cindx[m];
               g->ctype[m] = -1;
               m++;
               g->cindx[m] = -g->cindx[m-3];
               g->ctype[m] = -1;
               chainSize = m+1;
               m = 0;
               g->cindx[m] = -g->cindx[m+3];
               g->ctype[m] = -1;
               m = 1;
               g->cindx[m] = -g->cindx[m+1];
               g->ctype[m] = -1;
            }
         } //test=0

         if(g->test==1)
         {
            apply_periodic_LHS(&g[0],f1,f2,m); 
            chainSize = m+1;
            if(order==5) chainSize = m+2;
         }

      }

// ==================================================================

// ==================================================================
      for(j=0;j<chainSize;j++)
      {
         icell = g->cindx[j];
         itype = g->ctype[j];
         if (itype >=0) 
         {
            m = NVAR*icell;
            for(k=0;k<NVAR;k++)
            {
               consVar[k]=s->q[m];
               m++;
            }
            rhoi       = 1./consVar[0];
            g->f[j][0] = consVar[0];
            g->f[j][1] = consVar[1];
            g->f[j][2] = consVar[2];
            g->f[j][3] = consVar[3];
         }
         else
         {
            //
            // perform calculations on the ghost cells
            // based on whether they are on the solid boundary on that
            if (j < nghost) 
               iface = g->chainConn[f1];
            else
               iface = g->chainConn[f2-1];

            rightCell = g->faces[6*iface+4];
         
            // this is a face on solid wall 
            if (rightCell==-2)  
            {
               node1 = g->faces[6*iface];
               node2 = g->faces[6*iface+1];
        
               x1    = g->x[2*node1];
               y1    = g->x[2*node1+1];
               x2    = g->x[2*node2];
               y2    = g->x[2*node2+1];
        
               ds[0] = (y2-y1);
               ds[1] =-(x2-x1);

               icell =-icell;

               m = NVAR*icell;
               for(k=0;k<NVAR;k++)
               {
                  consVar[k] = s->q[m];
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
               g->f[j][0] = consVar[0];
               g->f[j][1] = (-consVar[1]*nx2ny-2*consVar[2]*nynx);
               g->f[j][2] = (consVar[2]*nx2ny-2*consVar[1]*nynx);
               g->f[j][3] = consVar[3];              
            }
            else 
            {
               g->f[j][0] = rinf;
               g->f[j][1] = rinf*s->uinf;
               g->f[j][2] = rinf*s->vinf;
               g->f[j][3] = pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);
            } // if loop (rightcell==-2)
         } // if loop (itype)
      } // j loop (chainsize)
      
      is  = nghost-1;
      ie  = chainSize-1;
      th  = 1./3.;
      qt  = 0.25;
      if (g->order==1) qt = 0.0;
      eps = 1e-10;

// ==================================================================
// Perform interface flux interpolation
// ==================================================================
      if(order==1 || order==3) 
      {
         muscld(g->f,g->ql,g->qr,g->f2,is,ie,th,qt,eps,chainSize,NVAR);
         // weno3(g->f,g->ql,g->qr,is,ie,eps,chainSize,NVAR);  //3rd weno
      }
      else if(order==5) 
         weno(g->f,g->ql,g->qr,is,ie,eps,chainSize,NVAR); 

    
      n   = is;
      idv = (g->chainConn[f1]==g->chainConn[f2-1]);

      for(f=f1;f<f2-idv;f++)
      {
         iface     = g->chainConn[f];
         node1     = g->faces[6*iface];
         node2     = g->faces[6*iface+1];
         
         leftCell  = g->faces[6*iface+2];
         rightCell = g->faces[6*iface+4];
         
         x1        = g->x[2*node1];
         y1        = g->x[2*node1+1];
         x2        = g->x[2*node2];
         y2        = g->x[2*node2+1];
         
         ds[0]     = (y2-y1);
         ds[1]     =-(x2-x1);

         for(m=0;m<NVAR;m++)
         {
            if (f==f2-idv-1 && idv==0) 
            {
               leftState[m]  = g->ql[n][m];
               rightState[m] = g->qr[n+1][m];
            }
            else
            {
               leftState[m]  = g->qr[n+1][m];
               rightState[m] = g->ql[n][m];
            }
         
            leftState0[m] = s->q[NVAR*leftCell+m]; //g->ql[j][m];               

            if (rightCell > -1 || g->test==1) 
               rightState0[m]=s->q[NVAR*rightCell+m]; //g->qr[j+1][m];
         }
     
         // periodic condition
         if (rightCell==-1 && g->test==0) 
         { 
            rightState0[0] = rinf;
            rightState0[1] = rinf*s->uinf;
            rightState0[2] = rinf*s->vinf;
            rightState0[3] = pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);
         }
         else if (rightCell==-2) 
         {
             dsnorm         = ds[0]*ds[0]+ds[1]*ds[1];
             nynx           = ds[0]*ds[1]/dsnorm;
             nx2ny          = (ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
             rightState0[0] = leftState0[0];
             rightState0[1] = -leftState0[1]*nx2ny-2*leftState0[2]*nynx;
             rightState0[2] = leftState0[2]*nx2ny-2*leftState0[1]*nynx;
             rightState0[3] = leftState0[3];
         }

// ==================================================================
// compute the Roe flux
// ==================================================================
     
         //
         //roeflx(&specRadius,flux,leftState,rightState,faceVel,ds,gm1);
         flux_roe2d_(ds,leftState,rightState,flux,&specRadius,&gamma1);
         jac_roe2d_(ds,leftState,rightState,lmat,rmat,&gamma1,&imode);

         //
         for(j=0;j<NVAR;j++)
            for(k=0;k<NVAR;k++)
            {
               (g->ff[iface]).lmat[j][k]=lmat[k][j];
               (g->ff[iface]).rmat[j][k]=rmat[k][j];
            }
         //
         m = NVAR*leftCell;
         for(j=0;j<NVAR;j++)
         {
            s->r[m]-=flux[j];
            m++;
         }
         s->sigma[leftCell]+=specRadius;
         if (rightCell > -1 || g->test==1) 
         {
            m=NVAR*rightCell;
            for(j=0;j<NVAR;j++)
            {
               s->r[m]+=flux[j];
               m++;
            }
            s->sigma[rightCell]+=specRadius;
         } 
         n++;
      } // f loop (f1 to f2)
   } // i loops (nchains)

   // compute L2 density residual
   *l2rho=0.;
   for(i=0;i<g->ncells;i++)
   {
      if ((*l2rho) < fabs(s->r[4*i])) 
      {
         icell=i;
         *l2rho=fabs(s->r[4*i]);
      }
   } // i loop 

}
// ##################################################################      
// END OF FILE
// ##################################################################
