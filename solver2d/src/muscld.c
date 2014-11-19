// ##################################################################
//
// muscld.c
//
// Koren's differentiable limiter for
// 3rd order accurate reconstruction
//
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include <stdlib.h>

void muscld(double **f,    // array of cons. quantities in the cell
            double **ql,   // left states
            double **qr,   // right states
            double **f2,   // temp flux states
            int      is,   // start index of chain
            int      ie,   // end index of chain
            double   th,   // third (=1/3)
            double   qt,   // quarter (=1/4)
            double   eps,  // epsilon in limiter
            int      imax, // chainSize
            int      nq)   // number of variables
{  
  int i,n,im;
  double thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt;
  
  im = imax-1;
  if (qt == 0.0)
  {
    for (n=0;n<nq;n++) // loop through number of variables
    {
      for (i=is;i<=ie;i++) // loop through the chain
      {
        ql[i][n] = f[i][n]; // left and right states
        qr[i][n] = f[i][n]; // are the same as the face
      }
    }
    return;
  } // if qt==0
  
  /*
   * do 3rd order otherwise
  */
  thm    = 1.0 - th;
  thp    = 1.0 + th;
 
  // loop through the states
  for(n=0;n<nq;n++) 
  {
    // compute f(i+1) - f(i)
    for(i=0;i<im;i++)
      f2[i][n]=f[i+1][n]-f[i][n];
 
    // loop from start to end of chain
    for(i=is;i<=ie;i++)
    {
      f2i      = f2[i][n];
      f2i1     = f2[i-1][n];
      a1       = 3.0*(f2i*f2i1+eps);
      a2       = 2.0*(f2i-f2i1)*(f2i-f2i1) + a1;
      f3i      = a1/a2 ;
      f3qt     = qt*f3i;
      ql[i][n] = f[i][n]+f3qt*( thm*f2i1 + thp*f2i );
      qr[i][n] = f[i][n]-f3qt*( thp*f2i1 + thm*f2i );
    } // i for loop
  } // n for loop
} // function end
// ##################################################################
// END OF FILE
// ##################################################################
