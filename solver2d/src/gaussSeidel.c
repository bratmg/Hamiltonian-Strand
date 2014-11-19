#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define NQ 4

void gaussSeidel(GRID *g,SOLN *s,double cflnum,double dt)
{
  //
  int i,j,k,l,m,f,n;
  int mm1,f1,f2,iface,idv,chainSize;
  double ds[2];
  double lmat[NQ][NQ];
  double rmat[NQ][NQ];
  double fmat[NQ];
  double dqmat[NQ];
  double Dmat[NQ][NQ];
  double dsnorm,nynx,nx2ny;
  int node1,node2,leftCell,rightCell,icell;
  double x1,y1,x2,y2,r01,r02,a1,a2,b1,b2,pp,dtfac;
  double linearl2rho;
  double damping;
  double tval[NQ];
  int isweep,jsweep;
  //
  // scale RHS and store
  // for linear iterations
  //
  for(i=0;i<g->ncells;i++)
    {
      
      dtfac=cflnum/s->sigma[i];
      for(m=0;m<NVAR;m++)
	{
         s->r[NVAR*i+m]*=dtfac;
	 s->r0[NVAR*i+m]=s->r[NVAR*i+m];
	 s->dq[NVAR*i+m]=0;
	}
    }
  //
  // find diagonal matrix terms
  //
  findDiagonals(g,s,cflnum,dt);
  //
  // now perform linear sweeps
  //
  for (isweep=0;isweep<g->msweep;isweep++) // number of linear sweeps
    {
      //
      // perform few Gauss-Seidel sweeps
      //  
      for(i=0;i<g->ncells*NVAR;i++)
	s->ddqf[i]=s->ddqb[i]=0;
      
      for (jsweep=0;jsweep<1;jsweep++) 
	{
	  for(i=0;i<g->ncells-1;i++)
	    {
	      for(j=0;j<NQ;j++)
		{
		  for(k=0;k<NQ;k++)
		    Dmat[j][k]=s->D[i][j][k];
		  fmat[j]=s->r[i*NVAR+j];
		}
	      //
	      invertMat4(Dmat,fmat,dqmat);
	      //
	      for(j=0;j<NQ;j++) s->ddq[NVAR*i+j]=dqmat[j];
	      //
	      for(j=0;j<4;j++)
		{
		  iface=g->c2f[4*i+j];
		  leftCell=g->faces[6*iface+2];
		  rightCell=g->faces[6*iface+4];
		  if (i==leftCell)
		    {
		      if (rightCell > -1 && rightCell > leftCell)
			{
			  dtfac=cflnum/s->sigma[rightCell];
			  axb((g->ff[iface]).lmat,&(s->ddq[NVAR*i]),&(s->ddqf[NVAR*i]),
			      &(s->r[NVAR*rightCell]),dtfac,NVAR);
			}
		    }
		  else if (leftCell > rightCell)
		    {
		      dtfac=cflnum/s->sigma[leftCell];
		      axb((g->ff[iface]).rmat,&(s->ddq[NVAR*i]),&(s->ddqf[NVAR*i]),
			  &(s->r[NVAR*leftCell]),-dtfac,NVAR);
		    }
		}
	    }
	  for(i=0;i<g->ncells*NVAR;i++)
	    s->ddqf[i]=s->ddq[i];
	  
	  for(i=g->ncells-1;i>=0;i--)
	    {
	      for(j=0;j<NQ;j++)
		{
		  for(k=0;k<NQ;k++)
		    Dmat[j][k]=s->D[i][j][k];
		  fmat[j]=s->r[i*NVAR+j];
		}
	      //
	      invertMat4(Dmat,fmat,dqmat);
	      //
	      for(j=0;j<NVAR;j++) s->ddq[i*NVAR+j]=dqmat[j];
	      //
	      for(j=0;j<4;j++)
		{
		  iface=g->c2f[4*i+j];
		  leftCell=g->faces[6*iface+2];
		  rightCell=g->faces[6*iface+4];
		  if (i==leftCell && rightCell < leftCell)
		    {
		      if (rightCell > -1)
			{
			  dtfac=cflnum/s->sigma[rightCell];
			  axb((g->ff[iface]).lmat,&(s->ddq[NVAR*i]),&(s->ddqb[NVAR*i]),
			      &(s->r[NVAR*rightCell]),dtfac,NVAR);
			}
		    }
		  else if (leftCell < rightCell)
		    {
		      dtfac=cflnum/s->sigma[leftCell];
		      axb((g->ff[iface]).rmat,&(s->ddq[NVAR*i]),&(s->ddqb[NVAR*i]),
			  &(s->r[NVAR*leftCell]),-dtfac,NVAR);
		    }
		}
	    }
	  for(i=0;i<g->ncells*NVAR;i++)
	    s->ddqb[i]=s->ddq[i];
	}
      
      for(i=0;i<g->ncells*NVAR;i++) s->dq[i]+=s->ddq[i];
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
      
      
  


