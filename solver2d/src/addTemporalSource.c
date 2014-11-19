#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void addTemporalSource(GRID *g, SOLN *s,double cflnum, double dt, int istep)
{
  int i,j,m;
  double dtau;
  double dtphys;

  if (g->timeInteg == BDF1 || istep==0) 
    {
      if (g->timeacc==0) 
	{
	  m=0;
	  for(i=0;i<g->ncells;i++)
	    {
             dtphys=cflnum*g->vol[i]/s->sigma[i];
	     for(j=0;j<NVAR;j++)
	       {
		 s->r[m]-=(s->q[m]-s->qt[m])*(g->vol[i]/dtphys);
		 m++;
	       }
	    }
	}
      if (g->timeacc==1) 
	{
	  m=0;
	  for(i=0;i<g->ncells;i++)
	    {
             dtphys=dt; 
	     for(j=0;j<NVAR;j++)
	       {
		 s->r[m]-=(s->q[m]-s->qt[m])*(g->vol[i]/dtphys);
		 m++;
	       }
	    }
	}
    }
  else if (g->timeInteg == BDF2 && istep > 0 ) 
    {
      if (g->timeacc==0) 
	{
	  m=0;
	  for(i=0;i<g->ncells;i++)
	    {
	      dtphys=TWOTHIRD*cflnum*g->vol[i]/s->sigma[i];
	      for(j=0;j<NVAR;j++)
		{
		  s->r[m]-=(s->q[m]-(2*s->qt[m]-0.5*s->qtt[m])*TWOTHIRD)*(g->vol[i]/dtphys);
		  m++;
		}
	    }
	}
      if (g->timeacc==1) 
	{
	  m=0;
	  for(i=0;i<g->ncells;i++)
	    {
	      dtphys=TWOTHIRD*dt;
	      for(j=0;j<NVAR;j++)
		{
		  s->r[m]-=(s->q[m]-(2*s->qt[m]-0.5*s->qtt[m])*TWOTHIRD)*(g->vol[i]/dtphys);
		  m++;
		}
	    }
	}
    }
}
 
