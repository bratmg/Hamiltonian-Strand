// ##################################################################
//
// stepSolution.c
//
// March through the solution for one time-step
// Choices of "step-type" are:
//  - euler, rk3, ADI, DADI, Gauss-Seidel, line-Gauss-Siedel
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <string.h>
#include <stdio.h>
#define deps 1e-10

void stepSolution(char *stepType,GRID *g,SOLN *s,double dt,double *l2norm)
{
  double     coef;
  static int istep=0;
  int        i,k,l;
  int        nsubiter=20;
  double     CFL;

  // create temp variable for CFL number
  CFL = g->CFL;

  if (strcmp(stepType,"euler")==0)
  {
    computeRHS(g,s,l2norm);
    coef = 1.0;
    updateSoln(s->q,s->q,s->r,s->sigma,dt,CFL,coef,g->ncells);
  }
  else if (strcmp(stepType,"rk3") == 0)
  {
    // RK step 1
    computeRHS(g,s,l2norm);
    coef=0.25;
    updateSoln(s->q,s->qt,s->r,s->sigma,dt,CFL,coef,g->ncells);
    coef=8.0/15;
    updateSoln(s->q,s->q,s->r,s->sigma,dt,CFL,coef,g->ncells);
    
    // RK step 2 
    computeRHS(g,s,l2norm);
    coef=5.0/12;
    updateSoln(s->qt,s->q,s->r,s->sigma,dt,CFL,coef,g->ncells);
    
    // RK step 3
    computeRHS(g,s,l2norm);
    coef=3.0/4.0;
    updateSoln(s->qt,s->q,s->r,s->sigma,dt,CFL,coef,g->ncells);
  }
  else if (strcmp(stepType,"ADI")==0)
  {
    computeRHSk(g,s,l2norm);
    ADI(g,s,s->cflnum,dt);
  }
  else if (strcmp(stepType,"DADI")==0)
  {
    computeRHSk(g,s,l2norm);
    DADI(g,s,s->cflnum,dt);
  }
  else if (strcmp(stepType,"Gauss-Seidel") == 0)
  {
    computeRHSkv(g,s,l2norm);
    gaussSeidel(g,s,s->cflnum,dt);
  }
  else if (strcmp(stepType,"line-Gauss-Seidel") == 0)
  {

    if (g->visc) 
      computeRHSkv(g,s,l2norm);
    else
      computeRHSk(g,s,l2norm);

    lineGaussSeidel(g,s,s->cflnum,dt);
  }
  istep++;
}
// ##################################################################
// END OF FILE
// ##################################################################
