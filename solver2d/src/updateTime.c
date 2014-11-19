#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void updateTime(GRID *g, SOLN *s)
{
  int i;

  for(i=0;i<NVAR*g->ncells;i++)
    {
      s->qtt[i]=s->qt[i];  // n -> n-1
      s->qt[i]=s->q[i];    // n+1 -> n
    }
}
 
