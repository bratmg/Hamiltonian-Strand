// ##############################################################
//
// weno.c
//
//
// Wtitten by Yong Su Jung
// ###############################################################
#include <stdio.h>
#include <stdlib.h>

void weno(double **f,  //fluxes at each edge
         double **ql, // left states
         double **qr, // right states
         int    is,   // start index of chain
         int    ie,   // end index of chain
         double epsw,  // epsilon in limiter
         int    imax, // chainSize
         int    nq)   // number of variables
{
 int i,n,im,k;
 double f0[imax-1],f1[imax-1],f2[imax-1];
 double slope[imax-1];
 double sol;



//  // 1st order upwind
//  for (n=0;n<nq;n++)
//  {
//   for (i=is;i<=ie;i++) //loop through chain 
//   {
//     ql[i][n] = f[i][n]; //left and right states are the same as the face
//     qr[i][n] = f[i][n];
//   }
//  }
//  return;
 
  // 5th order weno scheme
  im = imax-1;
  for (n=0;n<nq;n++) // number of flow variable
  {
   for (i=0;i<im;i++) // loop through the chain
   {   
     f0[i] = f[i][n];
   }
   for (i=is;i<=ie-1;i++) //1st difference at i+1/2
   {   
     f1[i] = f0[i+1] -f0[i];
   }
   for (i=is+1;i<=ie-1;i++) //2nd difference at i
   {
     f2[i] = f1[i] -f1[i-1];
   }


