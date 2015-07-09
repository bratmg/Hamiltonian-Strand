// ##################################################################
//
// driverMeshGen.c
//
// Main file that runs the mesh generation code
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meshtype.h"
#include "IO.h"
#include "globalVariables.h"
#include "listOperations.h"
#include "loopsAndCells.h"
#include "smoothOperations.h"
#include "memoryOperations.h"
#include "subdivision.h"

// ##################################################################
// begin main program
// ##################################################################
int main()
{

   int i,j;

   // welcome
   welcome();

	// read inputs
   readInputs();

   //
   // loop through grids
   //
   for (i = 0; i < nGrid; i++)
   {

   	// create output folder to store data in
   	createOutputFolder(i);

      // create list of triangles containing a given node
      createTriangleList(&g[i]);

      // find cell loops and vertex loops
      triCellAndVertexLoops(&g[i]);

      // colouring algorithm
      greedyColouringAlgorithm(&g[i]);

      // extract edges from triangulation 
      findTriangleEdges(&g[i]);
            
      // invert edge list (reverse list)
      createEdgeList(&g[i]);

      // create vertices on each edge (and quad edges)
      createVerticesOnEdge(&g[i]);

      // create quad cells
      createInteriorVertices(&g[i]);

   }

   //
   // loop through the hybrid grids and patch 
   // information to the regular grids
   //
   if(iHybrid)
   {
      for (i = 0; i < nHGrid; i++) meshPatching(&h[i]);
   }
   
   //
   // loop through grids
   //
   for (i = 0; i < nGrid; i++)
   {

      // find quad loops
      quadLoops(&g[i]);

      // smoothing
      if(strcmp(smoothTechnique,"lagrangian")==0)
      {
      	smoothGrid(&g[i],numSmooth);
		}
		else if(strcmp(smoothTechnique,"blend")==0)
		{
			smoothTriangleGrid(&g[i],numSmooth);
      	recreateVerticesOnEdge(&g[i]);
      	recreateInteriorVertices(&g[i]);
		}

      // create strand template
      if(iStrand) createStrands(&g[i]);

      // check quality and output statistics
      meshQuality(i,&g[i]);
      
      // write tecplot outputs
      writeTecplot(i,&g[i]);
   }

   //thanks
   thanks();

	return 0;	
}
// ##################################################################
// END OF FILE
// ##################################################################
