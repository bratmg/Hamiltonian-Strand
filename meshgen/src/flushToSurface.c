// ##################################################################
//
// flushToSurface.c
// 
// Routines that contain the routines to ensuee that the
// newly created points are along the surface of the boundary
// ##################################################################
#include <math.h>
#include <string.h>

#include "meshtype.h"
#include "flushToSurface.h"
#include "globalVariables.h"

#define INDX(row,col,ld) (row*ld + col)

void moveToBoundary(double *ptLeft, double *ptRight, 
                    double * pt, double *ptNormal, char *surfaceID)
{
   int    finiteTE;
   int    sign,s,cc,counter;
   double eps,invNormal,normal[2]; // only 2d
   double d,dNew,factor,delta;
   double t,c,x,y,z,x1,y1;
   double a0,a1,a2,a3,a4;

   double xNew, yNew, xtemp;

   // ===============================================================
   // For NACA airfoil
   // ===============================================================
   if(strcmp(surfaceID,"naca")==0 || strcmp(surfaceID,"NACA")==0)
   {
    

      x = pt[0];
      y = pt[1];

      t = 0.12;
      c = 1.00;
      x = x/c;

      a0 =  0.2969;
      a1 = -0.1260;
      a2 = -0.3516;
      a3 =  0.2843;

      finiteTE = 0;

      if(finiteTE == 1)
         a4 = -0.1015;
      else
         a4 = -0.1036;

      sign  = SIGN(pt[1]);

      // ==============================================================
      // Perform iterative technique to find the "normal" intersection
      // of the point with the airfoil
      // ==============================================================
      normal[0] = -(ptRight[1] - ptLeft[1]); // -(y2-y1)
      normal[1] =  (ptRight[0] - ptLeft[0]); //  (x2-x1)

      // ensure normal direction is accurate
      if(ptLeft[1]*normal[1] < 0)
      {
      // s = SIGN(ptLeft[1]*normal[1]);
         normal[0] = -normal[0];//(double)s*normal[0]; // -normal[0]
         normal[1] = -normal[1];//s*normal[1]; // -normal[1];
      }

      invNormal = 1./normal[1];

      eps = 1.;// initial 
      cc  = 0;
      // perform till convergence
      factor = 1.e-3;
      while (eps > 1.e-5)
      {
         d     = 0.;
         x1    = pt[0] + normal[0]*d;         
         xtemp = x1;
         while (x1 < 0. || x1 > 1.)
         {
            
            if(x1 > 1. && normal[0] > 0. ||
               x1 < 0. && normal[0] < 0.)
               x1 = (pt[0] - normal[0]*factor);   
            else 
               x1 = (pt[0] + normal[0]*factor);

            pt[0] = x1;
         }


         y1 = sign*5.0*t*c*( a0*sqrt(x1)   +
                             a1*x1         +
                             a2*x1*x1      +
                             a3*x1*x1*x1   +
                             a4*x1*x1*x1*x1  );

         // set new d
         dNew  = invNormal*(y1 - pt[1]);
         
         // error
         eps   = fabs(dNew);

         // update d
         d     = factor*dNew;

         // update position of point
         pt[0] = pt[0] + normal[0]*d;
         pt[1] = pt[1] + normal[1]*d;

         if(pt[0] != pt[0]) 
         {
            traces(nanning);
            printf("LeftPt:  (%e,%e)\n",ptLeft[0],ptLeft[1] );
            printf("RightPt: (%e,%e)\n",ptRight[0],ptRight[1] );
            printf("Normal:  (%e,%e)\n",normal[0],normal[1] );
            
            exit(1);
         }

      }
      
   }
   // ===============================================================
   // For a sphere
   // ===============================================================
   else if (strcmp(surfaceID,"sphere")==0)
   {
       // assumed radius of sphere as unity
       d = pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2];

       d = 1./sqrt(d);

       // renormalize the position
       pt[0] = pt[0]*d;
       pt[1] = pt[1]*d;
       pt[2] = pt[2]*d;
   }
   // ===============================================================
   // Robin Fuselage
   // ===============================================================
   else if (strcmp(surfaceID,"robin")==0)
   {
      double ptTemp[3];
      double sum;

      if (pt[0]>0.01 && pt[0]<1.97) 
      {
          robinSurface(pt);
          return;
      }


      ptTemp[0] = pt[0]; ptTemp[1] = pt[1]; ptTemp[2] = pt[2]; 

      eps = 1.;// initial 
      
      // perform till convergence
      factor = 0.001;
      counter = 0;
      delta = 1.e-6;
      d     = 0.;

      while (eps > 1.e-5 && counter < 20000)
      {
         
         ptTemp[0] = pt[0] + ptNormal[0]*d;
         ptTemp[1] = pt[1] + ptNormal[1]*d;
         ptTemp[2] = pt[2] + ptNormal[2]*d;


         if (ptTemp[0] < 0.)
         {
            ptTemp[0] = -ptTemp[0];
            delta*=0.1;
         }

         robinSurface(ptTemp);

         eps = fabs(ptTemp[1]-pt[1]);
      
         // update d
         d+=delta;


         if(pt[2] != pt[2]) 
         {
            traces(nanning);
            printf("LeftPt:  (%e,%e)\n",ptLeft[0],ptLeft[1] );
            printf("RightPt: (%e,%e)\n",ptRight[0],ptRight[1] );
            printf("Normal:  (%e,%e)\n",ptNormal[0],ptNormal[1] );
            
            exit(1);
         }
         counter++;
      }
      // printf("done.. eps (%e), counter(%d)\n",eps,counter );
      // traces(stopping in flush.);exit(1);
      pt[0] = ptTemp[0];
      pt[1] = ptTemp[1];
      pt[2] = ptTemp[2];
   }
   
}

// ==================================================================
//
// robinSurface
//
// Given a 'x' position, the function determines the 'y' and 'z'
// ==================================================================
void robinSurface(double *pt)
{
   int l = 8;
   int i,j,k;
   double phi,r,x,y,z,fs,fc;
   double hr,wr,zr,nr,temp,num,den,pwr;
   
   double cc[40] = { 1.0000, -1., -0.40, 0.40, 1.8, 0.00, 0.3825, 1.8000, 
                     0.3825,  0.,  0.00, 1.00, 0.0, 0.00, 1.0000, 1.0000, 
                     1.0000, -1., -0.96, 0.44, 3.0, 0.15, 0.2325, 0.1800,
                     0.1500,  0.,  0.00, 0.00, 0.0, 0.00, 1.0000, 1.0000, 
                     1.0000, -1., -1.9 , 0.1 , 2.0, 0.0 , 0.15  , 2.0000 };

   double dd[40] = { 1.0000, -1., -0.40, 0.40, 2.0, 0.00, 0.3275, 2.0000,
                     0.3275,  0.,  0.00, 1.00, 0.0, 0.00, 1.0000, 1.0000,
                     1.0000, -1., -0.96, 0.44, 3.0, 0.15, 0.1775, 0.1800,
                     0.1500,  0.,  0.00, 0.00, 0.0, 0.00, 1.0000, 1.0000,
                     1.0000, -1., -1.90, 0.10, 2.0, 0.00, 0.1500, 2.0000 };

   double ee[40] = { 1.0000, -1., -0.40, 0.40, 1.8, -0.0800,  0.0800, 1.8000,
                     0.0000,  0.,  0.00, 1.00, 0.0,  0.0000,  1.0000, 1.0000,
                     1.0000, -1., -0.96, 0.44, 3.0,  0.1163, -0.1163, 0.1800,
                     0.1163,  0.,  0.00, 0.00, 0.0,  0.0000,  1.0000, 1.0000,
                     0.1163,  0.,  0.00, 0.00, 0.0,  0.0000,  1.0000, 1.0000 };

   double ff[40] = { 2.,  3.,  0.00, 0.40, 1., 0., 1., 1.00,
                     5.,  0.,  0.00, 1.00, 0., 0., 1., 1.00,
                     1., -1., -0.96, 0.44, 1., 2., 3., 0.55,
                     2.,  0.,  0.00, 0.00, 0., 0., 1., 1.00,
                     2.,  0.,  0.00, 0.00, 0., 0., 1., 1.00 };

   rc = &cc[0];
   rd = &dd[0];
   re = &ee[0];
   rf = &ff[0];

   fs = 1.;
   fc = 1.;

   x = pt[0]; y = pt[1]; z = pt[2];

   // determine the coefficients to be used at each cross section
   k = 0;
   if (pt[0] >= 0.4 && pt[0] < 0.96)
      k = 1;
   else if (pt[0] >= 0.96 && pt[0] < 1.4)
      k = 2;
   else if (pt[0] >= 1.4 && pt[0] < 1.9)
      k = 3;
   else if (pt[0] >= 1.9 && pt[0] <= 2.0)
      k = 4;

   // obtain hr
   if (fabs(rc[INDX(k,3,l)]) < 1.e-8)
      temp = 1.;
   else
   {
      temp = fabs(pt[0]+rc[INDX(k,2,l)])/rc[INDX(k,3,l)];
      temp = pow(temp,rc[INDX(k,4,l)]);
   }

   temp = rc[INDX(k,0,l)] + rc[INDX(k,1,l)]*temp;
   temp = pow(temp,1./rc[INDX(k,7,l)]);

   hr = rc[INDX(k,5,l)] + rc[INDX(k,6,l)]*temp;


   // obtain wr
   if (fabs(rd[INDX(k,3,l)]) < 1.e-8)
      temp = 1.;
   else
   {
      temp = fabs(pt[0]+rd[INDX(k,2,l)])/rd[INDX(k,3,l)];
      temp = pow(temp,rd[INDX(k,4,l)]);
   }

   temp = rd[INDX(k,0,l)] + rd[INDX(k,1,l)]*temp;
   temp = pow(temp,1./rd[INDX(k,7,l)]);

   wr = rd[INDX(k,5,l)] + rd[INDX(k,6,l)]*temp;


   // obtain zr
   if (fabs(re[INDX(k,3,l)]) < 1.e-8)
      temp = 1.;
   else
   {
      temp = fabs(pt[0]+re[INDX(k,2,l)])/re[INDX(k,3,l)];
      temp = pow(temp,re[INDX(k,4,l)]);
   }
   
   temp = re[INDX(k,0,l)] + re[INDX(k,1,l)]*temp;
   temp = pow(temp,1./re[INDX(k,7,l)]);

   zr = re[INDX(k,5,l)] + re[INDX(k,6,l)]*temp;

   // obtain nr
   if (fabs(rf[INDX(k,3,l)]) < 1.e-8)
      temp = 1.;
   else
   {
      temp = fabs(pt[0]+rf[INDX(k,2,l)])/rf[INDX(k,3,l)];
      temp = pow(temp,rf[INDX(k,4,l)]);
   }

   temp = rf[INDX(k,0,l)] + rf[INDX(k,1,l)]*temp;
   temp = pow(temp,1./rf[INDX(k,7,l)]);

   nr = rf[INDX(k,5,l)] + rf[INDX(k,6,l)]*temp;

   // phi location of the point
   phi = atan(pt[1]/(pt[2]-zr));


   num = pow(0.25*hr*wr,nr);
   den = pow(0.5*hr*fabs(sin(phi)),nr) 
       + pow(0.5*wr*fabs(cos(phi)),nr);
   pwr = 1./nr;
   if(den > 0)
      r = pow(num/den,pwr);
   else
      r = 1.e-7;
   
   if(pt[2] < zr) {fs = -1.; fc = -1.;}

   pt[1] = r*fs*sin(phi);
   pt[2] = r*fc*cos(phi) + zr;

   // Test for NaN's
   if (pt[1] != pt[1] || pt[2] != pt[2])
   {
      printf("(x,y,z): (%e %e %e)\n",x,y,z);
      printf("(hr wr zr nr): (%e %e %e %e)\n",hr,wr,zr,nr);
      traces(Nanning in robin flush. Stopping.);
      exit(1);
   }



}

// ==================================================================
//
// createRobinData
//
// necessary arrays and variables for the robin fuselage
// ==================================================================
void createRobinData()
{


   // robin fuselage parameters
   // rc = (double *) malloc(sizeof(double)*40);
   // rd = (double *) malloc(sizeof(double)*40);
   // re = (double *) malloc(sizeof(double)*40);
   // rf = (double *) malloc(sizeof(double)*40);
   
   double cc[40] = { 1.0000, -1., -0.40, 0.40, 1.8, 0.00, 0.3825, 1.8000, 
                     0.3825,  0.,  0.00, 1.00, 0.0, 0.00, 1.0000, 1.0000, 
                     1.0000, -1., -0.96, 0.44, 3.0, 0.15, 0.2325, 0.1800,
                     0.1500,  0.,  0.00, 0.00, 0.0, 0.00, 1.0000, 1.0000, 
                     1.0000, -1., -1.9 , 0.1 , 2.0, 0.0 , 0.15  , 2.0000 };

   double dd[40] = { 1.0000, -1., -0.40, 0.40, 2.0, 0.00, 0.3275, 2.0000,
                     0.3275,  0.,  0.00, 1.00, 0.0, 0.00, 1.0000, 1.0000,
                     1.0000, -1., -0.96, 0.44, 3.0, 0.15, 0.1775, 0.1800,
                     0.1500,  0.,  0.00, 0.00, 0.0, 0.00, 1.0000, 1.0000,
                     1.0000, -1., -1.90, 0.10, 2.0, 0.00, 0.1500, 2.0000 };

   double ee[40] = { 1.0000, -1., -0.40, 0.40, 1.8, -0.0800,  0.0800, 1.8000,
                     0.0000,  0.,  0.00, 1.00, 0.0,  0.0000,  1.0000, 1.0000,
                     1.0000, -1., -0.96, 0.44, 3.0,  0.1163, -0.1163, 0.1800,
                     0.1163,  0.,  0.00, 0.00, 0.0,  0.0000,  1.0000, 1.0000,
                     0.1163,  0.,  0.00, 0.00, 0.0,  0.0000,  1.0000, 1.0000 };

   double ff[40] = { 2.,  3.,  0.00, 0.40, 1., 0., 1., 1.00,
                     5.,  0.,  0.00, 1.00, 0., 0., 1., 1.00,
                     1., -1., -0.96, 0.44, 1., 2., 3., 0.55,
                     2.,  0.,  0.00, 0.00, 0., 0., 1., 1.00,
                     2.,  0.,  0.00, 0.00, 0., 0., 1., 1.00 };

   rc = &cc[0];
   rd = &dd[0];
   re = &ee[0];
   rf = &ff[0];

}

// ##################################################################
// END OF FILE
// ##################################################################
