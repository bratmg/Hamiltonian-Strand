#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define NQ 5

void lineGaussSeidel1(GRID *g,SOLN *s,double cflnum,double dt)
{
  //
  int i,j,k,l,m,f,n,ii;
  int mm1;
  int f1,f2;
  int is,ie;
  int iface,iface1;
  int jface,jrightCell,jleftCell;
  int idv,idiff;
  int chainSize;
  int isweep;
  double damping;
  double ds[3];
  double leftState[NQ];
  double rightState[NQ];
  double leftState0[NQ];
  double rightState0[NQ];
  double consVar[NQ];
  double lmat[NQ][NQ];
  double rmat[NQ][NQ];
  double linearl2rho;
  double tval[NQ];
  double gm1=gamm-1.0;
  double specRadius;
  double linearL2rho;
  double l2rhodq;
  double faceVel=0.;
  double dsnorm,nynx,nx2ny;
  double rhoi;
  int node1,node2,leftCell,rightCell,icell;
  double x1,y1,x2,y2;  
  double r01,r02;
  double a1,a2,b1,b2;
  double pp;
  double th,qt,eps;
  double dscheck[2];
  double gamma1=gamm;
  double l2rho;
  int imode=1;
  double dtfac;
  FILE *fp1,*fp2,*fp3,*fp4,*fp5;
  double dtau;
  double *rdeltaq;
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  l2rho=0.;
  rdeltaq=(double *) malloc(sizeof(double)*NVAR);
  for(i=0;i<g->ncells;i++)
    {
      if ((l2rho) < fabs(s->r[4*i])) 
	  {
	    icell=i;
	    l2rho=fabs(s->r[4*i]);
	  }
      dtfac=cflnum/s->sigma[i];
      //dtfac=1;
      for(m=0;m<NVAR;m++)
	{
         s->r[NVAR*i+m]*=dtfac;
	 s->r0[NVAR*i+m]=s->r[NVAR*i+m];
	 s->dq[NVAR*i+m]=0;
	}
    }

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
      f1=g->faceStartPerChain[i];
      f2=g->faceStartPerChain[i+1];
      idv=(g->chainConn[f1]==g->chainConn[f2-1]);
      m=0;
      chainSize=(f2-idv-f1);
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
	  ds[0]=(y2-y1);
	  ds[1]=-(x2-x1);
	  ds[2]=0;
	  dsnorm=sqrt(ds[0]*ds[0]+ds[1]*ds[1]);
	  //
	  leftState0[0]=s->q[NVAR*leftCell];   
	  leftState0[1]=s->q[NVAR*leftCell+1]; 
	  leftState0[2]=s->q[NVAR*leftCell+2]; 
	  leftState0[3]=0;
	  leftState0[4]=s->q[NVAR*leftCell+3]; 
	  //
	  if (rightCell > -1) 
	    {
	      rightState0[0]=s->q[NVAR*rightCell]; 
	      rightState0[1]=s->q[NVAR*rightCell+1]; 
	      rightState0[2]=s->q[NVAR*rightCell+2]; 
	      rightState0[3]=0.;
	      rightState0[4]=s->q[NVAR*rightCell+3];       		      
	    }
	  //
	  if (rightCell==-1) {
	    rightState0[0]=rinf;
	    rightState0[1]=rinf*s->uinf;
	    rightState0[2]=rinf*s->vinf;
	    rightState0[3]=0;
	    rightState0[4]=pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);
	  }
	  else if (rightCell==-2) {
	    dsnorm=ds[0]*ds[0]+ds[1]*ds[1];
	    nynx=ds[0]*ds[1]/dsnorm;
	    nx2ny=(ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
	    rightState0[0]=leftState0[0];
	    rightState0[1]=-leftState0[1]*nx2ny-2*leftState0[2]*nynx;
	    rightState0[2]=leftState0[2]*nx2ny-2*leftState0[1]*nynx;
	    rightState0[3]=0;
	    rightState0[4]=leftState0[4];
	  }
	  //
	  jac_roe_(ds,leftState0,rightState0,lmat,rmat,&gamma1,&imode);	  
	  //
	  // transpose since jac_roe_ was in f90 ordering
	  // remove when jac_roe goes to C
	  //
	  for(j=0;j<NQ;j++)
	    for(k=0;k<NQ;k++)
             {
	      lmat[j][k]*=dsnorm;
              rmat[j][k]*=dsnorm;
             }

	  for(j=0;j<NQ;j++)
	    for(k=j+1;k<NQ;k++)
	      {
	  	swap(lmat[j][k],lmat[k][j]);
	  	swap(rmat[j][k],rmat[k][j]);
	     }
	  
	  if (rightCell > -1)
	    {
	      for(j=0;j<NQ;j++)
		for(k=0;k<NQ;k++)
		  {
		    s->D[leftCell][j][k]+=(lmat[j][k]);
		    s->D[rightCell][j][k]-=(rmat[j][k]);
		  }
	    }
	  else
	    {
	      if (rightCell==-2)
		{
		  dsnorm=ds[0]*ds[0]+ds[1]*ds[1];
		  nynx=ds[0]*ds[1]/dsnorm;
		  nx2ny=(ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
		  a1=-nx2ny;
		  b1=-2*nynx;
		  a2=-2*nynx;
		  b2=nx2ny;
		  for (j=0;j<NQ;j++)
		    {
		      r01=rmat[j][1];
		      r02=rmat[j][2];
		      rmat[j][1]=a1*r01+b1*r02;
		      rmat[j][2]=a2*r01+b2*r02;
		    }
		  for(j=0;j<NQ;j++)
		    for(k=0;k<NQ;k++)
		      s->D[leftCell][j][k]+=(rmat[j][k]);
		}
	      //
	      for(j=0;j<NQ;j++)
		for(k=0;k<NQ;k++)
		  s->D[leftCell][j][k]+=(lmat[j][k]);
	    }	  
	}
    }

  damping=1.0+0.5*cflnum*(g->msweep > 1);
  for(i=0;i<g->ncells;i++)
    {
      dtfac=cflnum/s->sigma[i];
      //dtfac=1./s->sigma[i];
      dtau = s->sigma[i]; //this is \frac{V}{\delta{\tau}}
      for(j=0;j<NQ;j++)
	for(k=0;k<NQ;k++)
	  {
	    s->D[i][j][k]*=dtfac;
	  }
      for(j=0;j<NQ;j++)
	s->D[i][j][j]+=damping;
      m++;
    }


  for (isweep=0;isweep<g->msweep;isweep++)
    {
      for (i=0;i<g->ncells*NVAR;i++) s->ddq[i]=0;

      for(ii=0;ii<2*g->nchains-1;ii++)
	{
	  i=(ii<g->nchains) ? ii: (2*g->nchains-2-ii);
	  idiff=(ii < g->nchains) ? 1: -1;
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
	  f1=g->faceStartPerChain[i];
	  f2=g->faceStartPerChain[i+1];
	  idv=(g->chainConn[f1]==g->chainConn[f2-1]);
	  m=0;
	  chainSize=(f2-idv-f1);
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
	      ds[0]=(y2-y1);
	      ds[1]=-(x2-x1);
	      ds[2]=0;
	      dsnorm=sqrt(ds[0]*ds[0]+ds[1]*ds[1]);
	      //
	      leftState0[0]=s->q[NVAR*leftCell];   
	      leftState0[1]=s->q[NVAR*leftCell+1]; 
	      leftState0[2]=s->q[NVAR*leftCell+2]; 
	      leftState0[3]=0;
	      leftState0[4]=s->q[NVAR*leftCell+3]; 
	      //
	      if (rightCell > -1) 
		{
		  rightState0[0]=s->q[NVAR*rightCell]; 
		  rightState0[1]=s->q[NVAR*rightCell+1]; 
		  rightState0[2]=s->q[NVAR*rightCell+2]; 
		  rightState0[3]=0.;
		  rightState0[4]=s->q[NVAR*rightCell+3];       		      
		}
	  
	      if (rightCell==-1) {
		rightState0[0]=rinf;
		rightState0[1]=rinf*s->uinf;
		rightState0[2]=rinf*s->vinf;
		rightState0[3]=0;
		rightState0[4]=pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);
	      }
	      else if (rightCell==-2) {
		dsnorm=ds[0]*ds[0]+ds[1]*ds[1];
		nynx=ds[0]*ds[1]/dsnorm;
		nx2ny=(ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
		rightState0[0]=leftState0[0];
		rightState0[1]=-leftState0[1]*nx2ny-2*leftState0[2]*nynx;
		rightState0[2]=leftState0[2]*nx2ny-2*leftState0[1]*nynx;
		rightState0[3]=0;
		rightState0[4]=leftState0[4];
	      }
	      //
	      jac_roe_(ds,leftState0,rightState0,lmat,rmat,&gamma1,&imode);	  
	      //
	      // transpose since jac_roe_ was in f90 ordering
	      // remove when jac_roe goes to C
	      //
	      for(j=0;j<NQ;j++)
		for(k=0;k<NQ;k++)
		  {
		    lmat[j][k]*=dsnorm;
		    rmat[j][k]*=dsnorm;
		  }

	      for(j=0;j<NQ;j++)
		for(k=j+1;k<NQ;k++)
		  {
		    swap(lmat[j][k],lmat[k][j]);
		    swap(rmat[j][k],rmat[k][j]);
		  }
	  
	      if (idv==1) {
		mm1=(m==0)?chainSize-1:m-1;
		for(j=0;j<NQ;j++)
		  for(k=0;k<NQ;k++)
		    {
		      g->A[m][j][k]+=(rmat[j][k]);
		      g->C[mm1][j][k]-=(lmat[j][k]);
		    }
		g->F[m][0]=s->r[NVAR*leftCell];
		g->F[m][1]=s->r[NVAR*leftCell+1];
		g->F[m][2]=s->r[NVAR*leftCell+2];
		g->F[m][3]=0.;
		g->F[m][4]=s->r[NVAR*leftCell+3];
	      }
	      else
		{
		  mm1=(m-1);
		  for(j=0;j<NQ;j++)
		    for(k=0;k<NQ;k++)
		      {
			if (mm1 > -1 && rightCell > -1)
			  {
			    g->A[m][j][k]+=(rmat[j][k]);
			    g->C[mm1][j][k]-=(lmat[j][k]);
			  }
		      }
		  //
		  if (rightCell < 0 && idv==0 && f==f2-1) m--;
		  //
		  g->F[m][0]=s->r[NVAR*leftCell];
		  g->F[m][1]=s->r[NVAR*leftCell+1];
		  g->F[m][2]=s->r[NVAR*leftCell+2];
		  g->F[m][3]=0.;
		  g->F[m][4]=s->r[NVAR*leftCell+3];
		}  
	      m++;
	    }     
	  m=0;
	  chainSize=chainSize-(idv==0);
	  for(f=f1;f<f2-1;f++)
	    {
	      iface=g->chainConn[f];
	      leftCell=g->faces[6*iface+2];
	      dtfac=cflnum/s->sigma[leftCell];
              //dtfac=1./s->sigma[leftCell];
	      //dtfac=0.;
	      for(j=0;j<NQ;j++)
		for(k=0;k<NQ;k++)
		  {
		    g->A[m][j][k]*=dtfac;
		    g->B[m][j][k]=s->D[leftCell][j][k];
		    g->C[m][j][k]*=dtfac;
		  }
	      //
	      iface1=g->chainConn[f+1];
	      icell=leftCell;
	      for(j=0;j<NVAR;j++) rdeltaq[j]=0;
	      //
	      for (j=0;j<4;j++)
		{
		  jface=g->c2f[4*icell+j];
		  if (jface != iface  &&
		      jface != iface1)
		    {
		      jleftCell=g->faces[6*jface+2];
		      jrightCell=g->faces[6*jface+4];
		      if (icell==jleftCell && jrightCell > -1) 
			{
			  axb1((g->ff[jface]).rmat,&(s->ddq[NVAR*jrightCell]),
			       rdeltaq,-dtfac,NVAR);
			}
		      else 
			{
			  axb1((g->ff[jface]).lmat,&(s->ddq[NVAR*jleftCell]),
			       rdeltaq,dtfac,NVAR);
			}
		    }
		}
	      g->F[m][0]+=rdeltaq[0];
	      g->F[m][1]+=rdeltaq[1];
	      g->F[m][2]+=rdeltaq[2];
	      g->F[m][3]+=0.;
	      g->F[m][4]+=rdeltaq[3];
	      m++;
	    }
	  //
	  // invert using appropriate banded block solver
	  //
	  if (idv==1) blockTridagPeriodic(g->A,g->B,g->C,g->F,chainSize,NQ);
	  if (idv==0) blockTridag(g->A,g->B,g->C,g->F,chainSize,NQ);
	  //
	  //
	  m=0;
	  for(f=f1;f<f2-1;f++)
	    {
	      //
	      // faces that bound the cell in 
	      // this chain
	      //
	      iface=g->chainConn[f];
	      iface1=g->chainConn[f+1];
	      //
	      icell=g->faces[6*iface+2];
	      //
	      s->ddq[NVAR*icell]  =g->F[m][0];
	      s->ddq[NVAR*icell+1]=g->F[m][1];
	      s->ddq[NVAR*icell+2]=g->F[m][2];
	      s->ddq[NVAR*icell+3]=g->F[m][4];
	      //
	      m++;
	    }
	}
      for(i=0;i<g->ncells*NVAR;i++) s->dq[i]+=s->ddq[i];
      l2rho=0.;
      icell=0;
      for(i=0;i<g->ncells;i++)
	{
	  if ((l2rho) < fabs(s->ddq[4*i])) 
	  {
	    icell=i;
	    l2rho=fabs(s->ddq[4*i]);
	  }
	}
      //trace(icell);
      //tracef(l2rho);
      computeLinearRHS(g,s,cflnum,&linearl2rho);
      tracef(linearl2rho);
    }

  l2rhodq=0.;
  for(i=0;i<g->ncells;i++)
    {
      if ((l2rhodq) < fabs(s->dq[4*i])) 
	{
	  icell=i;
	  l2rhodq=fabs(s->dq[4*i]);
	}
    }
  //trace(icell);
  tracef(l2rhodq);

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
  free(rdeltaq);
}
      
      
  


