// ##################################################################
//
// memoryOperations.h
//
// Memory related routines
// ##################################################################
#include <stdlib.h>

#include "meshtype.h"
#include "memoryOperations.h"

void freeMem(GRID *g)
{

   // free from loopsAndCells.c
   free(g->triLoop->index);
   free(g->triLoop->ID);
   free(g->triLoop);

   free(g->vertLoop->index);
   free(g->vertLoop->ID);
   free(g->vertLoop);





}
// ##################################################################
// END OF FILE
// ##################################################################