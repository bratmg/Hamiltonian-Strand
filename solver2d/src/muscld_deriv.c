// ##################################################################
// Koren's differentiable limiter for
// 3rd order accurate reconstruction
// ##################################################################
#include <stdio.h>
#include <stdlib.h>

void muscld_deriv(double **f,
                  double **ql,
                  double **qr,
                  double **f2,
                  double **dq,
                  int      is,
                  int      ie,
                  double   th,
                  double   qt,
                  double   eps,
                  int      imax,
                  int      nq)
{  
  int    i,n,im;
  double thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt;
  double df2i,df2i1,da1,da2,df3i,df3qt;

  im = imax-1;
  if (qt == 0.0)
  {
    for(n=0;n<nq;n++)
    {
      for(i=is;i<=ie;i++)
      {
        ql[i][n] =dq[i][n];
        qr[i][n] =dq[i][n];
      }
    }
    return;
  }
  //
  // do 3rd order otherwise
  //
  thm    = 1.0 - th;
  thp    = 1.0 + th;
  //
  // find the tangent
  //
  for(n=0;n<nq;n++)
  {
    for(i=0;i<im;i++) 
      f2[i][n]=f[i+1][n]-f[i][n];
    
    for(i=is;i<=ie;i++)
    {
      f2i      = f2[i][n];
      f2i1     = f2[i-1][n];
      df2i     = dq[i+1][n]-dq[i][n];
      df2i1    = dq[i][n]-dq[i-1][n];
      a1       = 3.0*(f2i*f2i1+eps);
      da1      = 3.0*(df2i*f2i1+f2i*df2i1);
      a2       = 2.0*(f2i-f2i1)*(f2i-f2i1) + a1;
      da2      = 2.0*2.0*(f2i-f2i1)*(df2i-df2i1)+da1;
      f3i      =  a1/a2 ;
      df3i     =  da1/a2 -a1*da2/(a2*a2);
      f3qt     = qt*f3i;
      df3qt    = qt*df3i;
      ql[i][n] = dq[i][n]+df3qt*( thm*f2i1 + thp*f2i ) + f3qt*(thm*df2i1+thp*df2i);
      qr[i][n] = dq[i][n]-df3qt*( thp*f2i1 + thm*f2i ) - f3qt*(thp*df2i1+thm*df2i);
    } // i loop
  } // n loop
}
// ##################################################################
// END OF FILE
// ##################################################################