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

void moveToBoundary(double *ptLeft, double *ptRight, 
                    double * pt, char *surfaceID)
{
   
   // ===============================================================
   // For NACA airfoil
   // ===============================================================
   if(strcmp(surfaceID,"naca")==0)
   {
      int    finiteTE;
      int    sign,cc;
      double eps,invNormal,normal[2]; // only 2d
      double d,dNew;
      double t,c,x,y,x1,y1;
      double a0,a1,a2,a3,a4;

      double xNew, yNew;


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
         normal[0] = -normal[0];
         normal[1] = -normal[1];
      }
      invNormal = 1./normal[1];

      eps = 1.;// initial 
      cc  = 0;
      // perform till convergence
      while (eps > 0.00001)
      {
         d  = 0.;
         x1 = pt[0] + normal[0]*d;         
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
         d     = 0.01*dNew;

         // update position of point
         pt[0] = pt[0] + normal[0]*d;
         pt[1] = pt[1] + normal[1]*d;

      }
      
   }
   // ===============================================================
   // For a sphere
   // ===============================================================
   else if (strcmp(surfaceID,"sphere")==0)
   {
       double d;

       // assumed radius of sphere as unity
       d = pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2];

       d = 1./sqrt(d);

       // renormalize the position
       pt[0] = pt[0]*d;
       pt[1] = pt[1]*d;
       pt[2] = pt[2]*d;
   }

   
}

// ##################################################################
// END OF FILE
// ##################################################################
