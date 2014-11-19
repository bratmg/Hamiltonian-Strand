#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include  "ham2dFunctionDefs.h"

void computeTimeScaling(GRID *g, SOLN *s,double cflnum,double dt,int istep)
{
  int i;
  double dtau;
  double dtphys;
  
  if (g->timeInteg==BDF1 || istep==0)
    {
      if (g-> timeacc== 0) 
	{
	  for(i=0;i<g->ncells;i++)
	    {
	      dtau=3.0*g->vol[i]/s->sigma[i]; // assume CFL = 1 for pseudo time step
	      dtphys=cflnum*dtau;         // physical time step
	      s->dtac[i]= (dtau*dtphys)/(dtau+dtphys)/g->vol[i]; // dt/dV
	      //s->dtac[i]= cflnum/s->sigma[i]; //(dtau*dtphys)/(dtau+dtphys)/g->vol[i]; // dt/dV
	    }
	}
      else if (g->timeacc==1)
	{
	  for(i=0;i<g->ncells;i++)
	    {
	      dtau=3.0*g->vol[i]/s->sigma[i]; // assume CFL = 1 for pseudo time step
	      dtphys=dt;                  // physical time step
	      s->dtac[i]= (dtau*dtphys)/(dtau+dtphys)/g->vol[i]; // dt/dV
              //s->dtac[i]=dtphys/g->vol[i];
	    }
	}
    }
  else if (g->timeInteg == BDF2 && istep > 0)
    {
      if (g-> timeacc== 0) 
	{
	  for(i=0;i<g->ncells;i++)
	    {
	      dtau=g->vol[i]/s->sigma[i]; // assume CFL = 1 for pseudo time step
	      dtphys=TWOTHIRD*cflnum*dtau;// physical time step with bdf2 scaling
	      s->dtac[i]= (dtau*dtphys)/(dtau+dtphys)/g->vol[i]; // dt/dV
	    }
	}
      else if (g->timeacc==1)
	{
	  for(i=0;i<g->ncells;i++)
	    {
	      dtau=g->vol[i]/s->sigma[i]; // assume CFL = 1 for pseudo time step
	      dtphys=TWOTHIRD*dt;         // physical time step with bdf2 scaling
	      s->dtac[i]= (dtau*dtphys)/(dtau+dtphys)/g->vol[i]; // dt/dV
	    }
	}
    }
}
 
