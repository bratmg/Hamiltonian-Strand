#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"


void computeRHSkv(GRID *g,SOLN *s,double *l2rho)
{
  //
  int i,j,k,m,f,n,kv;
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
  double leftStatev[NVAR];
  double rightStatev[NVAR];
  double leftState0[NVAR];
  double rightState0[NVAR];
  double lmat[NVAR][NVAR];
  double rmat[NVAR][NVAR];
  double lmatv[NVAR][NVAR];
  double rmatv[NVAR][NVAR];
  double consVar[NVAR];
  double flux[NVAR];
  double fluxv[NVAR];
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
  //
  // zero out residual and spectral radii
  //
  dscheck[0]=dscheck[1]=0;
  for(i=0;i<NVAR*g->ncells;i++) s->r[i]=0.0;
  for(i=0;i<g->ncells;i++) s->sigma[i]=0.0;
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  for(i=0;i<g->nchains;i++)
  {
    f1=g->faceStartPerChain[i];
    f2=g->faceStartPerChain[i+1];
    m=2;
    for(f=f1;f<f2;f++)
	{
	  iface=g->chainConn[f];
	  g->cindx[m]=g->faces[6*iface+2];
	  g->ctype[m]=1;
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
	  f=f1+1;
	  iface=g->chainConn[f];
	  g->cindx[m]=g->faces[6*iface+2];
	  g->ctype[m]=1;
	  chainSize=m;
	  m=0;
	  for(f=f2-3;f<f2-1;f++)
	    {
	      iface=g->chainConn[f];
	      g->cindx[m]=g->faces[6*iface+2];
	      g->ctype[m]=1;
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
	  m--;
	  g->cindx[m]=-g->cindx[m];
	  g->ctype[m]=-1;
	  m++;
	  g->cindx[m]=-g->cindx[m-2];
	  g->ctype[m]=-1;
	  chainSize=m+1;
	  m=0;
	  g->cindx[m]=-g->cindx[m+3];
	  g->ctype[m]=-1;
	  m=1;
	  g->cindx[m]=-g->cindx[m+1];
	  g->ctype[m]=-1;
	}
      for(j=0;j<chainSize;j++)
	{
	  icell=g->cindx[j];
	  itype=g->ctype[j];
	  if (itype >=0) 
	    {
	      m=NVAR*icell;
	      for(k=0;k<NVAR;k++)
		{
		  consVar[k]=s->q[m];
		  m++;
		}
	      rhoi=1./consVar[0];
	      //
	      // collect primitive variables in 
	      // g->f[j][0:3]
	      //
	      /*
	      g->f[j][0]=consVar[0];
	      g->f[j][1]=consVar[1]*rhoi;
	      g->f[j][2]=consVar[2]*rhoi;
	      g->f[j][3]=gm1*(consVar[3]-0.5*(consVar[1]*consVar[1]+consVar[2]*consVar[2])*rhoi); */
	      //if (i==610 && icell==19363) tracef(consVar[1]);
                
	      g->f[j][0]=consVar[0];
	      g->f[j][1]=consVar[1];
	      g->f[j][2]=consVar[2];
	      g->f[j][3]=consVar[3];	

	      g->fv[j][0]=consVar[0];
	      g->fv[j][1]=consVar[1];
	      g->fv[j][2]=consVar[2];
	      g->fv[j][3]=consVar[3];	

	    }
	  else
	    {
	      //
	      // do ghost cells
	      // based on whether they are on the solid boundary on that
	      if (j < 2) 
		{
		  iface=g->chainConn[f1];
		}
	      else
		{
		  iface=g->chainConn[f2-1];
		}
	      rightCell=g->faces[6*iface+4];
	      
	      if (rightCell==-2)  /* this is a face on solid wall */
		{
		  node1=g->faces[6*iface];
		  node2=g->faces[6*iface+1];
		  
		  x1=g->x[2*node1];
		  y1=g->x[2*node1+1];
		  x2=g->x[2*node2];
		  y2=g->x[2*node2+1];
		  
		  ds[0]=(y2-y1);
		  ds[1]=-(x2-x1);

		  icell=-icell;
		  m=NVAR*icell;
		  for(k=0;k<NVAR;k++)
		    {
		      consVar[k]=s->q[m];
		      m++;
		    }

		  dsnorm=ds[0]*ds[0]+ds[1]*ds[1];
		  nynx=ds[0]*ds[1]/dsnorm;
		  nx2ny=(ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
		  rhoi=1./consVar[0];
		  /*	
		  g->f[j][0]=consVar[0];
		  g->f[j][1]=(-consVar[1]*nx2ny-2*consVar[2]*nynx)*rhoi;
		  g->f[j][2]=(consVar[2]*nx2ny-2*consVar[1]*nynx)*rhoi;
		  g->f[j][3]=gm1*(consVar[3]-0.5*(consVar[1]*consVar[1]+consVar[2]*consVar[2])*rhoi);
		  */
		  g->f[j][0]=consVar[0];
		  g->f[j][1]=(-consVar[1]*nx2ny-2*consVar[2]*nynx);
		  g->f[j][2]=(consVar[2]*nx2ny-2*consVar[1]*nynx);
		  g->f[j][3]=consVar[3];		  		  
		  
		  g->fv[j][0]=consVar[0];
		  g->fv[j][1]=(-consVar[1]*nx2ny-2*consVar[2]*nynx);
		  g->fv[j][2]=(consVar[2]*nx2ny-2*consVar[1]*nynx);
		  g->fv[j][3]=consVar[3];		  		  


		}
	      else if (rightCell==-3) 
		{
		  icell=-icell;
		  m=NVAR*icell;
		  for(k=0;k<NVAR;k++)
		    {
		      consVar[k]=s->q[m];
		      m++;
		    }
		  
		  g->f[j][0]=consVar[0];
		  g->f[j][1]=0; //-consVar[1];
		  g->f[j][2]=0; //-consVar[2];		  
		  g->f[j][3]=consVar[3]
		    -(consVar[1]*consVar[1]+consVar[2]*consVar[2])*0.5/consVar[0];

		  g->fv[j][0]=consVar[0];
		  g->fv[j][1]=-consVar[1];
		  g->fv[j][2]=-consVar[2];		  
		  g->fv[j][3]=consVar[3];

		}
	      else
		{
		  g->f[j][0]=rinf;
		  g->f[j][1]=rinf*s->uinf;
		  g->f[j][2]=rinf*s->vinf;
		  g->f[j][3]=pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);

		  g->fv[j][0]=rinf;
		  g->fv[j][1]=rinf*s->uinf;
		  g->fv[j][2]=rinf*s->vinf;
		  g->fv[j][3]=pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);

		}
	    }
	}
      
      is=1;
      ie=chainSize-1;
      th=1./3;
      qt=0.25;
      if (g->order==1) qt=0.0;
      eps=1e-10;

      muscld(g->f,g->ql,g->qr,g->f2,is,ie,th,qt,eps,chainSize,NVAR);
//      weno3(g->f,g->ql,g->qr,is,ie,eps,chainSize,NVAR);  //3rd weno
//      weno(g->f,g->ql,g->qr,is,ie,eps,chainSize,NVAR);  //5th weno




      n=is;
      idv=(g->chainConn[f1]==g->chainConn[f2-1]);

      for(f=f1;f<f2-idv;f++)
	{
	  iface=g->chainConn[f];
	  node1=g->faces[6*iface];
	  node2=g->faces[6*iface+1];
	  leftCell=g->faces[6*iface+2];
	  rightCell=g->faces[6*iface+4];
	  x1=g->x[2*node1];
	  y1=g->x[2*node1+1];
	  x2=g->x[2*node2];
	  y2=g->x[2*node2+1];
	  //if (n==1) {
	  //  printf("%d %lf %lf %lf %lf\n",i,x1,x2,y1,y2);
	  //}
	  ds[0]=(y2-y1);
	  ds[1]=-(x2-x1);

	  for(m=0;m<NVAR;m++)
	    {
	      if (f==f2-idv-1 && idv==0) 
		{
		  leftState[m]=g->ql[n][m];
		  rightState[m]=g->qr[n+1][m];
		  leftStatev[m]=g->f[n][m];
		  rightStatev[m]=g->f[n+1][m];
		}
	      else
		{
		  leftState[m]=g->qr[n+1][m];
		  rightState[m]=g->ql[n][m];
		  leftStatev[m]=g->f[n+1][m];
		  rightStatev[m]=g->f[n][m];
		}
	      
	      leftState0[m]=s->q[NVAR*leftCell+m]; //g->ql[j][m];	            
	      if (rightCell > -1) 
		{
		  rightState0[m]=s->q[NVAR*rightCell+m]; //g->qr[j+1][m];
		}
	      
	    }
	  
          if (rightCell==-1) {
	    rightState0[0]=rinf;
	    rightState0[1]=rinf*s->uinf;
            rightState0[2]=rinf*s->vinf;
            rightState0[3]=pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);
          }
	  else if (rightCell==-2) {
	    dsnorm=ds[0]*ds[0]+ds[1]*ds[1];
	    nynx=ds[0]*ds[1]/dsnorm;
	    nx2ny=(ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
	    rightState0[0]=leftState0[0];
	    rightState0[1]=-leftState0[1]*nx2ny-2*leftState0[2]*nynx;
	    rightState0[2]=leftState0[2]*nx2ny-2*leftState0[1]*nynx;
	    rightState0[3]=leftState0[3];
	  }
	  
	  if (rightState[1]!=rightState0[1] && 0) {
	    trace(leftCell);
	    trace(rightCell);
            for(k=0;k<chainSize;k++)
              printf("%d %d %.16e\n",k,g->cindx[k],g->f[k][1]);
	    trace(n);
	    tracef(leftState[1]);
	    tracef(leftState0[1]);
	    tracef(rightState[1]);
	    tracef(rightState0[1]);
	    tracef(s->q[leftCell*NVAR+1])
	    trace(i);
            trace(n);
	    trace(f);
	    //exit(0);
	  }
	  //
	  if (rightCell < -1) 
	    {
	     wallflux_(ds,leftState,flux,&specRadius,&gamma1);
	     wallfluxjacobian_(ds,leftState,lmat,rmat,&gamma1);
	    }
	  else
	    {
	      flux_roe2d_(ds,leftState,rightState,flux,&specRadius,&gamma1);
	      jac_roe2d_(ds,leftState,rightState,lmat,rmat,&gamma1,&imode);
	   }
	  //
	  
	  //flux_visc_2d_(ds,&(g->vol[leftCell]),leftStatev,rightStatev,fluxv,
	  //	       &(s->rey),&(s->pr),&(s->prtr),&gamma1,&(s->c2b),
	  //	       &(s->rgas));
	  jac_visc_2d_(ds,&(g->vol[leftCell]),leftStatev,rightStatev,fluxv,lmatv,rmatv,
	  	       &(s->rey),&(s->pr),&(s->prtr),&gamma1,&(s->c2b),
	  	       &(s->rgas));
          
	  for(j=0;j<NVAR;j++)
            for(k=0;k<NVAR;k++)
              {
	        (g->ff[iface]).lmat[j][k]=(lmat[k][j]-lmatv[k][j]);
	        (g->ff[iface]).rmat[j][k]=(rmat[k][j]-rmatv[k][j]);
              }
	  //
	  m=NVAR*leftCell;
	  for(j=0;j<NVAR;j++)
	    {
	      s->r[m]-=flux[j];
	      s->r[m]+=fluxv[j];
	      m++;
	    }
	  s->sigma[leftCell]+=specRadius;
	  if (rightCell > -1) 
	    {
	      m=NVAR*rightCell;
	      for(j=0;j<NVAR;j++)
               {
		s->r[m]+=flux[j];
		s->r[m]-=fluxv[j];
                m++;
               }
	      s->sigma[rightCell]+=specRadius;
	    } 
	  n++;
	}
    }

  *l2rho=0.;
  for(i=0;i<g->ncells;i++)
    {
      *l2rho+=(s->r[4*i]*s->r[4*i]);
      //if ((*l2rho) < fabs(s->r[4*i])) 
	//  {
	 //   icell=i;
	  //  *l2rho=fabs(s->r[4*i]);
	  //}
    }
   *l2rho/=g->ncells;
   *l2rho=sqrt(*l2rho);
  //tracef(s->sigma[0]);
  //trace(icell);
}
      
      
  


