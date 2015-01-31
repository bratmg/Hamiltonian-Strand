// ##################################################################
//
// driverMeshGen.c
//
// Main file that runs the mesh generation code
// ##################################################################
#include <stdio.h>
#include <stdlib.h>

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

   int i;

   // welcome
   welcome();

	// read inputs
   readInputs();

   // loop through grids
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

      // find quad loops
      quadLoops(&g[i]);

      // smoothing
      if(0)
      {
      	smoothGrid(&g[i],numSmooth);
		}
		else
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
