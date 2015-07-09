// ##################################################################
//
// driverMeshGen.c
//
// Main file that runs the mesh generation code
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "IO.h"
#include "globalVariables.h"
#include "listOperations.h"
#include "loopsAndCells.h"
#include "smoothOperations.h"
#include "memoryOperations.h"
#include "subdivision.h"
#include "communication.h"

// ##################################################################
// begin main program
// ##################################################################
int main(int argc, char **argv)
{
   int i,ierr,nproc,myid;

   // start MPI communication
   ierr = MPI_Init(&argc,&argv);
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);

   // welcome
   if(myid == 0) welcome();
 
   // read inputs
   readInputs(myid);

   ierr = MPI_Barrier(MPI_COMM_WORLD);

   // create output folder to store data in
   createOutputFolder(myid);
   
   // create list of triangles containing a given node
   createTriangleList(&g[myid]);

   // find cell loops and vertex loops
   triCellAndVertexLoops(&g[myid]);

   // colouring algorithm
   greedyColouringAlgorithm(&g[myid]);

   // extract edges from triangulation 
   findTriangleEdges(&g[myid]);
            
   // invert edge list (reverse list)
   createEdgeList(&g[myid]);

   // create vertices on each edge (and quad edges)
   createVerticesOnEdge(&g[myid]);

   // create quad cells
   createInteriorVertices(&g[myid]);

   // find quad loops
   quadLoops(&g[myid]);

   // smoothing
   if(strcmp(smoothTechnique,"lagrangian")==0)
   {
   	smoothGrid(&g[myid],numSmooth);
	}
	else if(strcmp(smoothTechnique,"blend")==0)
	{
		smoothTriangleGrid(&g[myid],numSmooth);
     	recreateVerticesOnEdge(&g[myid]);
     	recreateInteriorVertices(&g[myid]);
	}

   // create boundary information
   boundaryNodeConnection_MPI(myid); 

   // create strand template
   if(iStrand) 
   {      
      // wait till thread 0 has finished
      ierr = MPI_Barrier(MPI_COMM_WORLD);

      createStrands(myid, &g[myid]);
   }

   // check quality and output statistics
   meshQuality(myid,&g[myid]);   
      
   // write tecplot outputs
   writeTecplot(myid,&g[myid]);

   ierr = MPI_Barrier(MPI_COMM_WORLD);

   //thanks
   if(myid == 0) thanks();

   MPI_Finalize();
	return 0;	
}
// ##################################################################
// END OF FILE
// ##################################################################
