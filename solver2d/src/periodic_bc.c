#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define nearBodyRadius 3.0

void periodic_bc(GRID *g)
{
  int i,m,i1,i2,i3,i4,kk;
  int f,f1,f2;
  int icell,iface;
  FILE *fp;
  int j,k,ff1,ff2,tmp;
  double bf1n[10000],bf2n[10000],bf3n[10000],bf4n[10000];

//////////////////////////////////////////////////////////////////
  // Find far-boundary face 
  //
  //       3
  //   |-------|
  //   |       |
  // 4 |       |2
  //   ---------
  //       1
  //  

  g->m1=0;
  g->m2=0;
  g->m3=0;
  g->m4=0;
  for(i=0;i<g->nfaces;i++)
  {
    if(g->faces[6*i+4] == -1)
    {
     i1 = g->faces[6*i];
     i2 = g->faces[6*i+1];
     if(g->x[2*i1+1]==g->x[2*i2+1]) // x parallel
     {
        if(g->x[2*i1+1] > 1.2) //3
        {
           g->bf3[g->m3] = i;
           bf3n[g->m3] = 0.5*(g->x[2*i1]+g->x[2*i2]);
            
           g->m3 = g->m3+1;
        }
        else //1
        {
           g->bf1[g->m1] = i;
           bf1n[g->m1] = 0.5*(g->x[2*i1]+g->x[2*i2]);
           g->m1 = g->m1+1;
        } 
     }
     if(g->x[2*i1] == g->x[2*i2]) // y parallel
     {
        if(g->x[2*i1]>1.2) //2
        {
           g->bf2[g->m2] = i;
           bf2n[g->m2] = 0.5*(g->x[2*i1+1]+g->x[2*i2+1]);
           g->m2 = g->m2+1;
        }
        else //4
        {
           g->bf4[g->m4] = i;
           bf4n[g->m4] = 0.5*(g->x[2*i1+1]+g->x[2*i2+1]);
           g->m4 = g->m4+1;
        }
     } //y parallel
    } //face = -1

  } // nface

     //////////////////////////////////////
     // re-locatinon in increasing order//
    //////////////////////////////////////

    // m=1
    for(i=0;i<g->m1;i++)
    {
     for(j=i+1;j<g->m1;j++)
     {
      if(bf1n[i]>bf1n[j])
      {
       tmp = g->bf1[i];
       g->bf1[i] = g->bf1[j];
       g->bf1[j] = tmp;
      }
     }
    }
    
    // m=2
    for(i=0;i<g->m2;i++)
    {
     for(j=i+1;j<g->m2;j++)
     {
      if(bf2n[i]>bf2n[j])
      {
       tmp = g->bf2[i];
       g->bf2[i] = g->bf2[j];
       g->bf2[j] = tmp;
      }
     }
    }

    // m=3
    for(i=0;i<g->m3;i++)
    {
     for(j=i+1;j<g->m3;j++)
     {
      if(bf3n[i]>bf3n[j])
      {
       tmp = g->bf3[i];
       g->bf3[i] = g->bf3[j];
       g->bf3[j] = tmp;
      }
     }
    }

    // m=4
    for(i=0;i<g->m4;i++)
    {
     for(j=i+1;j<g->m4;j++)
     {
      if(bf4n[i]>bf4n[j])
      {
       tmp = g->bf4[i];
       g->bf4[i] = g->bf4[j];
       g->bf4[j] = tmp;
      }
     }
    }

   /////////////////////////////////////
   //make a neighbors face information// 
   /////////////////////////////////////
    // m1
    for (j=0;j<g->m1;j++)
    {
    for(k=0;k<g->nchains;k++)
    {
       f1 = g->faceStartPerChain[k];
       f2 = g->faceStartPerChain[k+1];
       ff1=g->chainConn[f1];
       ff2=g->chainConn[f2-1];
       if(g->bf1[j] == ff1) 
       {
        g->bf1_neigf[j][0] = g->chainConn[f1+1]; 
        g->bf1_neigf[j][1] = g->chainConn[f1+2]; 
        g->bf1_neigf[j][2] = g->chainConn[f1+3]; 
        break;
       } 
       if(g->bf1[j] ==ff2) 
       {
        g->bf1_neigf[j][0] = g->chainConn[f2-2]; 
        g->bf1_neigf[j][1] = g->chainConn[f2-3]; 
        g->bf1_neigf[j][2] = g->chainConn[f2-4]; 
        break;
       }
    }
//    printf("face index = %d, n_face_index = %d %d %d\n",g->bf1[j],g->bf1_neigf[j][0],g->bf1_neigf[j][1],g->bf1_neigf[j][2]);
    }

    // m2
    for (j=0;j<g->m2;j++)
    {
    for(k=0;k<g->nchains;k++)
    {
       f1 = g->faceStartPerChain[k];
       f2 = g->faceStartPerChain[k+1];
       ff1=g->chainConn[f1];
       ff2=g->chainConn[f2-1];
       if(g->bf2[j] == ff1) 
       {
        g->bf2_neigf[j][0] = g->chainConn[f1+1]; 
        g->bf2_neigf[j][1] = g->chainConn[f1+2]; 
        g->bf2_neigf[j][2] = g->chainConn[f1+3]; 
        break;
       } 
       if(g->bf2[j] ==ff2) 
       {
        g->bf2_neigf[j][0] = g->chainConn[f2-2]; 
        g->bf2_neigf[j][1] = g->chainConn[f2-3]; 
        g->bf2_neigf[j][2] = g->chainConn[f2-4]; 
        break;
       } 
    }
    }

    // m3
    for (j=0;j<g->m3;j++)
    {
    for(k=0;k<g->nchains;k++)
    {
       f1 = g->faceStartPerChain[k];
       f2 = g->faceStartPerChain[k+1];
       ff1=g->chainConn[f1];
       ff2=g->chainConn[f2-1];
       if(g->bf3[j] == ff1) 
       {
        g->bf3_neigf[j][0] = g->chainConn[f1+1]; 
        g->bf3_neigf[j][1] = g->chainConn[f1+2]; 
        g->bf3_neigf[j][2] = g->chainConn[f1+3]; 
        break;
       } 
       if(g->bf3[j] ==ff2) 
       {
        g->bf3_neigf[j][0] = g->chainConn[f2-2]; 
        g->bf3_neigf[j][1] = g->chainConn[f2-3]; 
        g->bf3_neigf[j][2] = g->chainConn[f2-4]; 
        break;
       } 
    }
    }

    // m4
    for (j=0;j<g->m4;j++)
    {
    for(k=0;k<g->nchains;k++)
    {
       f1 = g->faceStartPerChain[k];
       f2 = g->faceStartPerChain[k+1];
       ff1=g->chainConn[f1];
       ff2=g->chainConn[f2-1];
       if(g->bf4[j] == ff1) 
       {
        g->bf4_neigf[j][0] = g->chainConn[f1+1]; 
        g->bf4_neigf[j][1] = g->chainConn[f1+2]; 
        g->bf4_neigf[j][2] = g->chainConn[f1+3]; 
        break;
       } 
       if(g->bf4[j] ==ff2) 
       {
        g->bf4_neigf[j][0] = g->chainConn[f2-2]; 
        g->bf4_neigf[j][1] = g->chainConn[f2-3]; 
        g->bf4_neigf[j][2] = g->chainConn[f2-4]; 
        break;
       } 
    }

    }
}
