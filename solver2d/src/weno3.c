// ##############################################################
//
// weno3.c
//
//
// Wtitten by Yong Su Jung
// ###############################################################
#include <stdio.h>
#include <stdlib.h>

void weno3(double **f, //fluxes at each edge
         double **ql, // left states
         double **qr, // right states
         int    is,   // start index of chain
         int    ie,   // end index of chain
         double epsw,  // epsilon in limiter
         int    imax, // chainSize
         int    nq)   // number of variables

{
  int i,n,im,k;
  double f2[imax-1][nq];
  double f2i,f2i1,a1,a2,w01,w02,w11,w12,f3;
 
  im = imax-1; 

  for (n=0;n<nq;n++)  //number of flow variable
  {
   for (i=0;i<im;i++) //loop through the chain
   {   
     f2[i][n] = f[i+1][n]-f[i][n]; // difference
   }

   for (i=is;i<im;i++) //calcualte
   {   
     f2i      = f2[i][n];
     f2i1     = f2[i-1][n];
     a1       = 1.0/(f2i1*f2i1+epsw);
     a2       = 1.0/(f2i*f2i+epsw);
     w01      = a1/(a1+2.0*a2);
     w02      = 1 - w01;
     w11      = 2.0*a1/(2.0*a1+a2);
     w12      = 1 -w11;
     ql[i][n] = f[i][n]+0.5*( w01*f2i1 + w02*f2i);
     qr[i][n] = f[i][n]-0.5*( w11*f2i1 + w12*f2i);

   }
   //first-order at boundary
     qr[0][n] = f[0][n];
     ql[0][n] = f[0][n];

     qr[im][n] = f[im][n];
     ql[im][n] = f[im][n];

   
   }




} //function end

   

 











        


