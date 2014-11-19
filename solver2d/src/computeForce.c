// ##################################################################
//
// computeForce.c
//
// ##################################################################
#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void computeForce(GRID *g,SOLN *s)
{
  int i,n;
  double x1,x2,y1,y2;
  int iface,node1,node2,icell;
  double rho,rhou,rhov,e,p;
  double fac;
  double cs,ss;

  s->cx=s->cy=s->cl=s->cd=0;

  for (i=0;i<g->nbfaces;i++)
  {
    iface  = g->bfaces[i];
    node1  = g->faces[6*iface];
    node2  = g->faces[6*iface+1];
    
    x1     = g->x[2*node1];
    y1     = g->x[2*node1+1];
    x2     = g->x[2*node2];
    y2     = g->x[2*node2+1];
    
    icell  = g->faces[6*iface+2];
    //
    rho    = s->q[NVAR*icell];
    rhou   = s->q[NVAR*icell+1];
    rhov   = s->q[NVAR*icell+2];
    e      = s->q[NVAR*icell+3];
    //
    p      = (gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov)/rho);
    //
    s->cy -= p*(x2-x1);
    s->cx += p*(y2-y1);
  } // i loop

  fac    = 0.5*rinf*s->mach*s->mach;
  s->cy /= fac;
  s->cx /= fac;

  cs     = cos(s->alpha*deg2rad);
  ss     = sin(s->alpha*deg2rad);
  s->cl  = s->cy*cs-s->cx*ss;
  s->cd  = s->cy*ss+s->cx*cs;      
}
// ##################################################################
// END OF FILE
// ##################################################################