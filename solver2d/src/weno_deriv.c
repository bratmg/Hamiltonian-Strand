/**
   WENO's differentiable
*/

#include <stdio.h>
#include <stdlib.h>

double weno5_d(double a, double b, double c, double d, double e, double a_d, double b_d, double c_d, double d_d, double e_d,double epsw);


void weno_deriv(double **f,
		double **ql,
		double **qr,
		double **dq,
		int is,
		int ie,
		double eps,
		int imax,
		int nq)
{  
  int i,n,im;
  double epsf;
  double f0[imax],df0[imax];

  im = imax-2;
  for (n=0;n<nq;n++)
  {
   for(i=0;i<=ie;i++) // loop through the chain
   {
     f0[i] = f[i][n];
     df0[i] = dq[i][n];
   }

   for (i=is;i<im;i++)
   {
     ql[i][n]=weno5_d(f0[i-2],f0[i-1],f0[i],f0[i+1],f0[i+2],df0[i-2],df0[i-1],df0[i],df0[i+1],df0[i+2],eps);
     qr[i][n]=weno5_d(f0[i+2],f0[i+1],f0[i],f0[i-1],f0[i-2],df0[i+2],df0[i+1],df0[i],df0[i-1],df0[i-2],eps);
   }
   }
}



//##########################################################
//
// weno5_d
//
//##########################################################
//#include <stdio.h>
//#include <stdlib.h>

double weno5_d(double a,  //fluxes at each edge
          double  b, // left states
          double  c, // right states
          double  d,   // start index of chain
          double  e,   // end index of chain
          double  a_d,
          double  b_d, // left states
          double  c_d, // right states
          double  d_d,   // start index of chain
          double  e_d,   // end index of chain
          double  epsw)
{
 double b1,b2,djm1,ejm1,dj,ej,djp1,ejp1;
 double dis0,dis1,dis2,q30,q31,q32,d01,d02,a1ba0,a2ba0;
 double w0,w1,w2;
 double sol;
 
 double djm1_d,ejm1_d,dj_d,ej_d,djp1_d,ejp1_d;
 double dis0_d,dis1_d,dis2_d,q30_d,q31_d,q32_d,d01_d,d02_d,a1ba0_d,a2ba0_d;
 double w0_d,w1_d,w2_d;

      
      b1 = 13./12.;
      b2 = 1./6.;
      djm1 = a-2.*b+c;
      djm1_d = a_d-2.*b_d+c_d;


      ejm1 = a-4.*b+3.*c;
      ejm1_d = a_d-4.*b_d+3.*c_d;

      dj   = b-2.*c+d;
      dj_d = b_d-2.*c_d+d_d;

      ej   = b-d;
      ej_d = b_d - d_d;

      djp1 = c-2.*d+e;
      djp1_d = c_d -2.*d_d+e_d;


      ejp1 = 3.*c-4.*d+e;
      ejp1_d = 3.*c_d-4.*d_d+e_d;

      dis0 = b1*djm1*djm1+0.25*ejm1*ejm1+epsw;
      dis0_d = b1*2.*djm1*djm1_d+0.25*2.*ejm1_d;

      dis1 = b1*dj*dj+0.25*ej*ej+epsw;
      dis1_d = b1*2.*dj*dj_d + 0.25*2.*ej*ej_d;

      dis2 = b1*djp1*djp1+0.25*ejp1*ejp1+epsw;
      dis2_d = b1*2.*djp1*djp1_d + 0.25*2.*ejp1*ejp1_d;

      q30 = 2.*a-7.*b+11.*c;
      q30_d = 2.*a_d - 7.*b_d + 11.*c_d;

      q31 = -b+5.*c+2.*d;
      q31_d = -b_d+5.*c_d+2.*d_d;

      q32 = 2.*c+5.*d-e;
      q32_d = 2.*c_d + 5.*d_d -e_d;

      d01 = dis0/dis1;
      d01_d = (dis0_d*dis1-dis0*dis1_d)/(dis1*dis1);

      d02 = dis0/dis2;
      d02_d = (dis0_d*dis2-dis0*dis2_d)/(dis2*dis2);

      a1ba0 = 6.*d01;
      a1ba0_d = 6.*d01_d;

      a2ba0 = 3.*d02;
      a2ba0_d = 3.*d02_d;

      w0 = 1.0/(1.0+a1ba0+a2ba0);
      w0_d = -(a1ba0_d+a2ba0_d)/((1.0+a1ba0+a2ba0)*(1.0+a1ba0+a2ba0));

      w1 = a1ba0*w0;
      w1_d = a1ba0_d*w0 + a1ba0*w0_d;

      w2 = 1.-w0-w1;
      w2_d = -w0_d-w1_d;

//      sol = b2*(w0*q30+w1*q31+w2*q32);  //output

      sol = b2*(w0_d*q30 + w0*q30_d + w1_d*q31+w1*q31_d + w2_d*q32 + w2*q32_d);

      return sol;
        
        
} // end function
