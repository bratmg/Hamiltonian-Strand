#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
void apply_periodic(GRID *g,int f1, int f2, int m)
{  
  int i,j,k,f,n;
  int is,ie;
  int iface;
  int idv;
  double leftState[NVAR];
  double rightState[NVAR];
  double leftState0[NVAR];
  double rightState0[NVAR];
  int node1,node2,node3,node4,leftCell,rightCell,icell;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;  
  double xa,ya,za,xb,yb,zb;
  double pp;
  double th,qt,eps;
  int iflag;
  int order, nghost;
  
  order = g->order;
               nghost = 2;
  if(order==5) nghost = 3;

  iface = g->chainConn[f1];

  for(j=0;j<g->m1;j++)
  {
    if(iface == g->bf1[j]) 
    {g->cindx[nghost-1] = g->faces[6*g->bf3[j]+2];
     g->cindx[nghost-2] = g->faces[6*g->bf3_neigf[j][1]+2];
     if(order==5) 
     {
      g->cindx[nghost-3] = g->faces[6*g->bf3_neigf[j][2]+2]; 
     }
     break;
    }
    
    if(iface == g->bf3[j]) 
    {g->cindx[nghost-1] = g->faces[6*g->bf1[j]+2];
     g->cindx[nghost-2] = g->faces[6*g->bf1_neigf[j][1]+2];            
     if(order==5) 
     {
      g->cindx[nghost-3] = g->faces[6*g->bf1_neigf[j][2]+2]; 
     }            
     break;
    }
   }


  for(j=0;j<g->m2;j++)
  {
    if(iface == g->bf2[j]) 
    {g->cindx[nghost-1] = g->faces[6*g->bf4[j]+2];
     g->cindx[nghost-2] = g->faces[6*g->bf4_neigf[j][1]+2];
     if(order==5) 
     {
      g->cindx[nghost-3] = g->faces[6*g->bf4_neigf[j][2]+2]; 
     }            
     break;
    }
    if(iface == g->bf4[j]) 
    {g->cindx[nghost-1] = g->faces[6*g->bf2[j]+2];
     g->cindx[nghost-2] = g->faces[6*g->bf2_neigf[j][1]+2];
     if(order==5) 
     {
      g->cindx[nghost-3] = g->faces[6*g->bf2_neigf[j][2]+2]; 
     }           
     break;
    }
  }


  iface = g->chainConn[f2-1];
  for(j=0;j<g->m1;j++)
  {
    if(iface == g->bf1[j]) 
    {g->cindx[m-1] = g->faces[6*g->bf3[j]+2];
     g->cindx[m] = g->faces[6*g->bf3_neigf[j][0]+2];
     if(order==5) 
     {
      g->cindx[m+1] = g->faces[6*g->bf3_neigf[j][1]+2]; 
     }
     break;
    }
    if(iface == g->bf3[j]) 
    {g->cindx[m-1] = g->faces[6*g->bf1[j]+2];
     g->cindx[m] = g->faces[6*g->bf1_neigf[j][0]+2];
     if(order==5) 
     {
      g->cindx[m+1] = g->faces[6*g->bf1_neigf[j][1]+2]; 
     }
     break;
    }
  }           

  for(j=0;j<g->m2;j++)
  {
    if(iface == g->bf2[j]) 
    {g->cindx[m-1] = g->faces[6*g->bf4[j]+2];
     g->cindx[m] = g->faces[6*g->bf4_neigf[j][0]+2];
     if(order==5) 
     {
      g->cindx[m+1] = g->faces[6*g->bf4_neigf[j][1]+2]; 
     }
     break;
    }
    if(iface == g->bf4[j]) 
    {g->cindx[m-1] = g->faces[6*g->bf2[j]+2];
     g->cindx[m] = g->faces[6*g->bf2_neigf[j][0]+2];
     if(order==5) 
     {
      g->cindx[m+1] = g->faces[6*g->bf2_neigf[j][1]+2]; 
     }
     break;
    }
  }


}
//END of SUBROUTINE  
