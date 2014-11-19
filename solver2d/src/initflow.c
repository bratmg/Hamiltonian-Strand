// ##################################################################
//
// initflow.c
//
// Initialize the flow 
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NQ 4

void initflow(GRID *g,SOLN *s)
{
  int i,j,m;
  FILE *fp;

  double xc,yc;
  int n1,n2,n3,n4;
  double x1,y1,x2,y2,x3,y3,x4,y4;
  double vstr,xv0,yv0;
  double pi,c1,c2,xx,yy,r2,rr;
  double r,p,du,dv,tp,dtp,tp0;
  double gm1;
  double xv,yv;

  printf("init_test:%d\n",g->test);

  //
  fp = fopen("input.ham2d","r");
  fscanf(fp,"Mach=%lf\n",&s->mach);
  fscanf(fp,"alpha=%lf\n",&s->alpha);
  fscanf(fp,"rey=%lf\n",&s->rey);
  fclose(fp);

  s->rey=s->rey/s->mach; // Mach scaled Reynolds number

  tracef(s->mach);
  tracef(s->alpha);
  
  // Mach scaled velocities?
  s->uinf = s->mach*cos(s->alpha*deg2rad);
  s->vinf = s->mach*sin(s->alpha*deg2rad);
  s->einf = pinf/(gamm-1)+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);

  tracef(s->uinf);
  tracef(s->vinf);
  
  s->gm1  = gamm-1.;
  s->c2b  = 0.3678;
  s->rgas = 1./gamm;
  s->pr   = 0.72;
  s->prtr = 0.3333;

  //
  // allocate all the solution arrays
  //
  s->q    = (double *)   malloc(sizeof(double)*NVAR*g->ncells);
  s->qt   = (double *)   malloc(sizeof(double)*NVAR*g->ncells);
  s->qtt  = (double *)   malloc(sizeof(double)*NVAR*g->ncells);
  s->r    = (double *)   malloc(sizeof(double)*NVAR*g->ncells);
  s->dq   = (double *)   malloc(sizeof(double)*NVAR*g->ncells);
  s->ddq  = (double *)   malloc(sizeof(double)*NVAR*g->ncells);
  s->ddqb = (double *)   malloc(sizeof(double)*NVAR*g->ncells);
  s->ddqf = (double *)   malloc(sizeof(double)*NVAR*g->ncells);
  s->dtac = (double *)   malloc(sizeof(double)*g->ncells);
  s->r0   = (double *)   malloc(sizeof(double)*NVAR*g->ncells);
  s->D    = (double ***) malloc(sizeof(double **)*g->ncells);
  s->itag = (int *)      malloc(sizeof(int)*g->ncells);
  g->ff   = (faceMat *)  malloc(sizeof(faceMat)*g->nfaces);
  
  //
  for (i=0;i<g->ncells;i++)
  {
    s->D[i]=(double **) malloc(sizeof(double *)*NQ);
      
    for (j=0;j<NQ;j++)
      s->D[i][j]=(double *) malloc(sizeof(double)*NQ);
  }
  
  s->sigma=(double *) malloc(sizeof(double)*g->ncells);

  // Loop through the number of cells and initialize 
  // freestream conservative variables for q and qt
  m=0;

// ==================================================================
// far-field bc
// ==================================================================  
  if(g->test==0)
  {
    for (i=0;i<g->ncells;i++)
    {
      s->q[m]  = rinf;
      s->qt[m] = s->q[m];
      m++;

      s->q[m]  = rinf*s->uinf;
      s->qt[m] = s->q[m];
      m++;

      s->q[m]  = rinf*s->vinf;
      s->qt[m] = s->q[m];
      m++;
      
      // why can't the next line be = s->einf?
      s->q[m]  = pinf/(gamm-1.) + 0.5*rinf*(s->uinf*s->uinf+
                  s->vinf*s->vinf);
      s->qt[m] = s->q[m];
      m++;
    } // i loop
  }

// ==================================================================  
// for isentropic vortex problem
// ==================================================================  
  if(g->test == 1) {
    for (i=0;i<g->ncells;i++)
    {
      n1 = g->conn[4*i];
      n2 = g->conn[4*i+1];
      n3 = g->conn[4*i+2];
      n4 = g->conn[4*i+3];

      xc = (g->x[2*n1]+g->x[2*n2]+g->x[2*n3]+g->x[2*n4])/4.;
      yc = (g->x[2*n1+1]+g->x[2*n2+1]+g->x[2*n3+1]+g->x[2*n4+1])/4.;
      
      pi  = deg2rad*180.;
      gm1 = gamm-1.;
      tp0 = 1.0/gamm;

      s->uinf = 0.5;
      s->vinf = 0.0;
      vstr    = 0.0;//5.0;
      xv0     = 5.;
      yv0     = 5.;
   
   
      c1  = 0.5*vstr/pi;
      c2  = 0.125*gm1/(gamm*gamm)*pow(vstr,2)/pow(pi,2);//gamm^2?
    
      xv  = 0.;
      yv  = 0.;
    
      xx  = xc - (xv0+xv);
      yy  = yc - (yv0+yv);
    
      r2  = xx*xx + yy*yy;
      rr  = sqrt(r2);
    
      du  = c1*exp(0.5*(1.-r2))*(-yy);
      dv  = c1*exp(0.5*(1.-r2))*(xx);
      dtp = -c2*exp(1.-r2);
    
      tp = tp0 + dtp;
    
      r = pow(gamm*tp,1.0/gm1);
   
      s->uinf  = s->uinf + du;
      s->vinf  = s->vinf + dv;
      p        = pow(r,gamm);
    
      s->q[m]  = r;
      s->qt[m] = s->q[m];
      m++;

      s->q[m]  = r*s->uinf;
      s->qt[m] = s->q[m];
      m++;

      s->q[m]  = r*s->vinf;
      s->qt[m] = s->q[m];
      m++;

      s->q[m]  = p/(gamm-1.)+0.5*r*(s->uinf*s->uinf+s->vinf*s->vinf); 
      s->qt[m] = s->q[m]; //what is this?

      m++;    
    } // i loop
  } // test==1

}
// ##################################################################  
// END OF FILE
// ##################################################################
