// ##################################################################
//
// outputSolution.c
//
// Subroutines for output plots
//
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

// ==================================================================
// 
// ==================================================================
void outputSolution(GRID *g,SOLN *s,int nn)
{
  int    i,n;
  int    iface,node1,node2;
  int    icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE   *fp,*fp1;
  char fname[80];

  //

//  fp=fopen("./output/output.plt","w");
  if      (nn < 10  ) {sprintf(fname,"./output/output00%d.dat",nn);}
  else if (nn < 100 ) {sprintf(fname,"./output/output0%d.dat",nn);}
  else if (nn < 1000) {sprintf(fname,"./output/output%d.dat",nn);}
  fp=fopen(fname,"w");

  //fp1=fopen("check.dat","w");
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\",\"RES\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED)\n");

  // list of x-coordinate positions (after smoothing)
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");

  // list of y-coordinate positions (after smoothing)
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");

  // rho, rhoU, rhoV, E for all the cells
  for(n=0;n<NVAR;n++)
    for (i=0;i<g->ncells;i++)
    {
       fprintf(fp,"%f\n",s->q[4*i+n]);
    }
  fprintf(fp,"\n");
  // Added by Bharath (Testing. Residual.)
  for (i=0;i<g->ncells;i++)
  {
    fprintf(fp,"%f\n",fabs(s->r[4*i]));
  }    
  // End added by Bharath (Testing. Residual.)

  

  // connectivity information for the cells
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);

  fclose(fp);
  fp=fopen("./output/cp.dat","w");
  for(i=0;i<g->nbfaces;i++)
   {
    iface = g->bfaces[i];
    node1 = g->faces[6*iface];
    node2 = g->faces[6*iface+1];
    x1    = g->x[2*node1];
    y1    = g->x[2*node1+1];
    x2    = g->x[2*node2];
    y2    = g->x[2*node2+1];
    icell = g->faces[6*iface+2];
    rho   = s->q[NVAR*icell];
    rhou  = s->q[NVAR*icell+1];
    rhov  = s->q[NVAR*icell+2];
    e     = s->q[NVAR*icell+3];
    pp    = (gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov)/rho);
    cp    = (pp-pinf)/(0.5*s->mach*s->mach);  

    fprintf(fp,"%.16e %.16e %.16e\n",(x1+x2)*0.5,(y1+y2)*0.5,cp);
  }
  fclose(fp);
  //fclose(fp1);
}

// ==================================================================
// 
// ==================================================================
void outputdq(GRID *g,SOLN *s)
{
  int i,n;
  int iface,node1,node2;
  int icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE *fp;
  char fname[80];
  static int istep0=0;
  //
  sprintf(fname,"dq%d.plt",istep0);
  fp=fopen(fname,"w");
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->ncells;i++)
      fprintf(fp,"%f\n",s->ddq[4*i+n]);
  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);
  fclose(fp);
  istep0++;
}

// ==================================================================
// 
// ==================================================================
void outputr(GRID *g,SOLN *s)
{
  int i,n;
  int iface,node1,node2;
  int icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE *fp;
  char fname[80];
  static int istep0=0;
  //
  sprintf(fname,"r%d.plt",istep0);
  fp=fopen(fname,"w");
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->ncells;i++)
      fprintf(fp,"%f\n",s->r[4*i+n]);
  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);
  fclose(fp);
  istep0++;
}
// ##################################################################
// END OF FILE
// ##################################################################s
