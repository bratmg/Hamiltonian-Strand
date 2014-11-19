// ##################################################################
//
// ham2d.c
//
// Driver code for the 2d unstrctured Hamiltonian path based code
//
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <time.h>
#include <math.h>

int main()
{
  GRID    *g;     // pointer to GRID data structure
  SOLN    *s;     // pointer to SOLN data structure
  int     ngrids; // total number of grids
  int     nsteps; // total number of time steps
  int     nwrite; // total number of time steps
  int     i,n,nn;    //
  double  dt;     // time step
  double  CFL;
  double  l2rho;
  int     msweep;
  char    c;
  char    fname[20];
  char    scheme[20];
  char    timeInteg[20];
  int     order,timeacc,visc,test;
  FILE    *fp;
  clock_t start, end;
  double  cpu_time_used;
  
// ==================================================================  
// Read in file inputs from input.ham2d
// ==================================================================
  fp=fopen("input.ham2d","r");
  while((c=fgetc(fp))!='\n'); // skip the 
  while((c=fgetc(fp))!='\n'); // first three
  while((c=fgetc(fp))!='\n'); // lines of input.ham2d
  fscanf(fp,"scheme=%s\n",scheme);
  fscanf(fp,"time integration=%s\n",timeInteg);
  fscanf(fp,"order=%d\n",&order);
  fscanf(fp,"timeacc=%d\n",&timeacc);
  fscanf(fp,"nsteps=%d\n",&nsteps);
  fscanf(fp,"nwrite=%d\n",&nwrite);
  fscanf(fp,"dt=%lf\n",&dt);
  fscanf(fp,"CFL=%lf\n",&CFL);
  fscanf(fp,"msweep=%d\n",&msweep);
  fscanf(fp,"visc=%d\n",&visc);
  fscanf(fp,"testcase=%d\n",&test);
  fclose(fp);
  trace(nsteps);
  tracef(dt);
  //
  ngrids=1;
  g=(GRID *) malloc(sizeof(GRID)*ngrids);
  s=(SOLN *) malloc(sizeof(SOLN)*ngrids);

// ==================================================================
// preprocess grids
// code is written with an overset
// framework in mind (but not implemented yet)
// (ngrid==1) for now
// ==================================================================
  for(i=0;i<ngrids;i++) 
  {
      // Read the grid data (i.e., data points and connectivity information
      // from the Matlab mesh-generation)
      readGrid(&g[i]);
      g[i].visc = visc;

      // Run the preprocessing steps (computes cell vol., and
      // obtains cell to face and cell to chain connectivity)
      g[i].test = test;
      preprocess(&g[i]);

      // Initialize the flow
      initflow(&g[i],&s[i]);
      // Initialize other grid and solution variables
      g[i].order   = order;   // order of the scheme
      g[i].timeacc = timeacc; // time accuracy
      g[i].CFL     = CFL;     // CFL number for grid
      s[i].cflnum  = CFL;     // CFL number for soln
      g[i].msweep  = msweep;  // number of sweeps 

      if (strcmp(timeInteg,"bdf1")==0) g[i].timeInteg=BDF1;
      if (strcmp(timeInteg,"bdf2")==0) g[i].timeInteg=BDF2;
    }

  tracef(CFL)

// ==================================================================
// Main bulk of the code. Step through the solution for 
// as many steps are required.
// ==================================================================
  fp = fopen("./output/sol_his.dat","w");
  
  printf("#ham2d : using %s scheme for inversion\n",scheme);
  printf("#ham2d : using %d order of spatial accurate for invicid\n",order);
  fprintf(fp,"ham2d : using %s scheme for inversion\n",scheme);
  fprintf(fp,"ham2d : using %d order of spatial accurate for invicid\n",order);
  
  cpu_time_used=0.;

  for (n=0;n<nsteps;n++) // loop through the time steps
  { 
    for(i=0;i<ngrids;i++) // loop for all the grids
    { 
      // start the clock
      start = clock();

      // Step through the solution for one time step
      stepSolution(scheme,&g[i],&s[i],dt,&l2rho);

      // compute the force on the airfoil
      computeForce(&g[i],&s[i]);

      // end the clock
      end = clock();

      // compute CPU time used
      cpu_time_used += (((double) (end - start)) / CLOCKS_PER_SEC);
      
      printf("%d %e %2.4f %2.4f %2.4f\n",n,l2rho,s[i].cl,s[i].cd,cpu_time_used);
      fprintf(fp,"%d %e %2.4f %2.4f %2.4f\n",n,l2rho,s[i].cl,s[i].cd,cpu_time_used);

    } // i loop

    // Output the solution to screen and files
    if((n+1)%nwrite==0 || n==0) 
    {
      nn = (n+1)/nwrite;
      outputSolution(&g[0],&s[0],nn); 
    }
  } // n loop

   fclose(fp);

}
// ##################################################################
// END OF FILE
// ##################################################################
