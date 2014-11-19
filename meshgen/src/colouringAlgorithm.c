// ##################################################################
//
// colouringAlgorithm.c
// 
// Contains the various colouring algorithms
// ##################################################################
#include <stdio.h>
#include <stdlib.h>

#include "meshtype.h"

#define COLOUR_CHAR_MAX 64
#define NC 7

// ##################################################################
// colour the loops
// ##################################################################
void greedyColouringAlgorithm(GRID *g)
{

   printf("#meshgen: Greedy colouring algorithm...\n");

   int i,j,k,kk;
   int lenLoop;
   int loop[g->vertLoop->maxLen];
   int okcol[g->vertLoop->maxLen];

   g->colourChar  = (char *) malloc(sizeof(char)*COLOUR_CHAR_MAX);
   g->colourID    = (int *)  malloc(sizeof(int)*g->numTriNode);
   g->colourIndex = (int *)  malloc(sizeof(int)*g->numTriNode);

   // initialize the colour arrays
   i = 0;
   while (i+NC < COLOUR_CHAR_MAX)
   {
      g->colourChar[i  ] = 'b';
      g->colourChar[i+1] = 'r';
      g->colourChar[i+2] = 'g';
      g->colourChar[i+3] = 'm';
      g->colourChar[i+4] = 'c';
      g->colourChar[i+5] = 'y';
      g->colourChar[i+6] = 'k';
      i += NC;
   }

   for (i = 0; i < g->numTriNode; i++) {g->colourID[i] = -1;}

   for (i = 0; i < g->vertLoop->maxLen; i++) {okcol[i] = 10000;}

   g->maxCol = 0;

   // loop over total number of triangular nodes
   for (i = 0; i < g->numTriNode; i++)
   {
      // vertex loop centered around node i
      k = 0;
      for (j = g->vertLoop->index[i]; j < g->vertLoop->index[i+1]; j++)
      {         
         loop[k] = g->vertLoop->ID[j];
         k++;
      }

      lenLoop = g->vertLoop->index[i+1] - g->vertLoop->index[i];
      
      // okcol = (1:maxcol)
      for (j = 0; j <= g->maxCol; j++) { okcol[j] = j; }

      for (k = 0; k < lenLoop; k++)
      {           
         if (g->colourID[loop[k]] > -1 )
         {            
            okcol[g->colourID[loop[k]]] = g->maxCol + 1;
         }
      } // j loop

      // colr(i) = min(okcol)      
      g->colourID[i] = findMinMax(okcol,g->vertLoop->maxLen,'min');
      g->maxCol      = MAX(g->maxCol,g->colourID[i]);


   } // i loop

   // initialize numCol
   g->numCol = (int *) malloc(sizeof(int)*(g->maxCol+1));

   kk = 0;
   // loop through all the nodes
   for (j = 0; j <= g->maxCol; j++)
   {
      k = 0;   
      for (i = 0; i < g->numTriNode; i++)      
      {   
         if( g->colourID[i] == j) 
         {
            g->colourIndex[kk] = i;
            kk++; // counter for colourIndex
            k++;  // counter for number of nodes per colour
         }
      }
      
      // number of nodes of a particular colour
      g->numCol[j] = k;      
   }

 //  writearrayINT(g->numCol,g->numTriNode);
}


// ##################################################################
// find min or max of an int array
// ##################################################################
int findMinMax(const int *array, const int size, const char *option)
{
   int i;
   int min,max;
   min =  10000;
   max = -10000;
   
   if (option == 'min')
   {
      for (i = 0; i < size; i++){min = MIN(min,array[i]);}

      return min;
   }
   else if (option == 'max')
   {
      for (i = 0; i < size; i++){max = MAX(max,array[i]);}

      return max;
   }
   else
   {
      traces(In colouringAlgorithm.c);
      traces(Incorrect option. Stopping.);
      exit(1);
   }

}


// ##################################################################
// END OF FILE
// ##################################################################
