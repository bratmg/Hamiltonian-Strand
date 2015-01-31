// ##################################################################
//
// listOperations.c
//
// ##################################################################
#include <stdio.h>
#include <stdlib.h>

#include "meshtype.h"
#include "IO.h"

// ##################################################################
//
// reverseList.c
//
//  Creates reverse list and associated index list
//
// list1 is a         list of elements of Set 2 per element of Set 1 (size N)
// list2 is a reverse list of elements of Set 1 per element of Set 2 (size N)
// n1 is the number of elements in Set 1
// n2 is the number of elements in Set 2
// index1 is the starting index in list1 for each element of Set 1 (size n1+1)
// index2 is the starting index in list2 for each element of Set 2 (size n2+1)
//
// ==================================================================
// E.g. list2 = [218,315,316||13,234,216||13,38,96||..]
//     index2 = [ 1, 4, 7, ... ] (need not be offset by 3)
//
// The lines in _conn.dat file with vertex ID i can be found in 
// list2 by looking at indiced from index(i) to index(i+1)-1
//
// Therefore, Vertex ID 2 is contained in lines 13, 246, and 216
// ==================================================================
//
// ##################################################################
void reverseList(int *list1,
                 int *list2,
                 int  N,
                 int  n1,
                 int  n2,
                 int *index1,
                 int *index2)
{

   int i,i2,j;
   int k1,k2[n2];
   int Ninst[n2];

   // n2 is the number of points in coord file
   for (i = 0; i < n2; i++)
      Ninst[i] = 0;


   // N is (3*nt) - total vertex list
   // the following loop builds a list of the number of
   // times a given vertex appears in the _conn.dat file.
   // Implies the number of edges connecting a given vertex
   // e.g.: if Ninst(230) = 7, then vertex ID 229 is connected to 7 edges
   for (i = 0; i < N; i++)
      Ninst[list1[i]]++;   

   // First point in index2 is set to 1
   index2[0]=0;

   // Loop from 2 to n2+1
   // n2 is the number of data points in coord file
   // index2[i] contains the index2[i-1] + the number of 
   // edges connected to vertex[i-1]
   for (i = 1; i < n2+1; i++)
      index2[i] = index2[i-1] + Ninst[i-1];   

   
   for (i = 0; i < n2; i++)
      k2[i] = -1;

   // Loop over n1 (lines in conn file)
   k1 = -1;
   for (i = 0; i < n1; i++)
   {
      
      // this loop runs from 1:3 or x:x+2
      for (j = index1[i]; j < index1[i+1]; j++)
      {      
         k1++; 

         // effectively goes thourgh all the elements in conn file
         // n1*3 ( 3 = size of j loop)
         i2 = list1[k1];
         k2[i2]++;

         list2[index2[i2] + k2[i2] ] = i;
         
      } // j loop
      
   } // iloop
   
}

// ##################################################################
//
// createTriangleList
//
// Given the starting and ending index of a triangle, the nodes of that
// triangle are known. The following function, however, outputs data
// such that, the triangles containing a given node are now listed
// ##################################################################
void createTriangleList(GRID *g)
{

   printf("#meshgen: Creating triangle list ...\n");

   int i;
   int node2triIndex[g->numTriangle+1];

   // allocate space for the list and index arrays
   g->tri2nodeList  = (int *) malloc(sizeof(int)*3*g->numTriangle);
   g->tri2nodeIndex = (int *) malloc(sizeof(int)*(g->numTriNode+1));
   

   node2triIndex[0] = 0;
   for (i = 1; i < g->numTriangle+1 ; i++)
      node2triIndex[i] = node2triIndex[i-1] + 3;

   // reverse list
   reverseList(g->triConn,g->tri2nodeList,                    // list1,list2
               3*g->numTriangle,g->numTriangle,g->numTriNode, // N,n1,n2
               node2triIndex,g->tri2nodeIndex);               // index1,index2


}

// ##################################################################
//
// createEdgeList
//
// Given an edge, using the starting and ending index, triangles the
// edges belong to are known. The following function, however, outputs data
// such that, the edge IDs of each triangle are now listed
// ##################################################################
void createEdgeList(GRID *g)
{

   printf("#meshgen: Creating edge list ...\n");

   int i,k;
   int *index1;
   int *list1, *listTemp;

   index1   = (int *) malloc(sizeof(int)*(g->numTriEdge+1));

   listTemp = (int *) malloc(sizeof(int)*(2*g->numTriEdge));

   k = 0; // running counter for list1
   index1[0] = 0;

   // loop over all identified triangle edges
   for (i = 0; i < g->numTriEdge; i++)
   {
      index1[i+1] = index1[i];

      // Edges belong to some triangle
      // (not sure when this will fail)
      if (g->triEdge[i][2] > -1)
      {
         listTemp[k] = g->triEdge[i][2];
         index1[i+1]++;
         k++;
      }

      // if edge shared by two triangles
      if (g->triEdge[i][3] > -1)
      {
         listTemp[k] = g->triEdge[i][3];
         index1[i+1]++;
         k++;
      }

   }
   // allocation and initialization of list1
   list1 = (int *) malloc(sizeof(int)*k);   
   for (i = 0; i < k; i++)
      list1[i] = listTemp[i];

   // allocate space for the list and index arrays
   g->edge2triList  = (int *) malloc(sizeof(int)*k);
   g->edge2triIndex = (int *) malloc(sizeof(int)*(g->numTriangle+1));
   

   // reverse list
   reverseList(list1,g->edge2triList,          // list1,list2
               k,g->numTriEdge,g->numTriangle, // N,n1,n2
               index1,g->edge2triIndex);       // index1,index2

}
// ##################################################################
// END OF FILE
// ##################################################################
