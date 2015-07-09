// ##################################################################
//
// meshtype.h
// 
// Contains the structure definition used in the code
// ##################################################################
#ifndef MESHTYPE_H
#define MESHTYPE_H

// #include <process.h>
// #include <dir.h>
#include <mpi.h>

// ==================================================================
typedef struct LOOP
{
   int  maxLen; // maximum length of any single loop
   int  totLen; // total length of all loops
   int *index;  // identified start and end indices of loops 
   int *ID;     // list of vertices/cells (their ID) in the loop

} LOOP;

// ==================================================================
typedef struct EDGE
{
   int   numEdge; // total number of triangular edges 
   int  *ID;      // 4/6 element format for each edge

} EDGE;

// ==================================================================
typedef struct STRAND
{   
   double  normal[3]; // normal direction of the strand
   double *pos;       // [x,y,z] position along a strand

} STRAND;


// ==================================================================
typedef struct GRID
{   
   
   // triangular grid
   int      numTriNode;    // original triangular nodes
   int      numTriangle;   // total number of triangles

   int     *triConn;       // connectivity data for the triangles
   double  *nodePosTri;    // coordinates of the nodes [x,y,z. ordered data]
   double  *nodeNormal; // outward normal vectors of the nodes

   int     *tri2nodeList;  // contains the list of triangles for every node
   int     *tri2nodeIndex; // contains the indexing of the list

   // various loops
	LOOP    *triLoop;       // loop for triangles around a node
	LOOP    *vertLoop;      // loop for the vertices around a node
   LOOP    *quadLoop;      // loop for the quadrilateral loops around a node
   LOOP    *strandLoop;    // "loop" for the various strands
   LOOP    *cellLoop;      // loop across all cells (every layer and strand)

   // colouring algorithm
   int      maxCol;        // maximum number of colours
   char    *colourChar;    // character array for the colours
   int     *colourID;      // contains the colour ID form colourChar for each node
   int     *colourIndex;   // node ID associated with the colour (mapped to colourID)
   int     *numCol;        // number of nodes of a particular colour

   // triangle edges
   int     numTriEdge;     // number of triangle edges
	int   **triEdge;        // possible array for triangle edge

   int     *edge2triList;  // contains the list of triangles for every edge
   int     *edge2triIndex; // contains the indexing of the list

   // quad cells and edges
   int      numQuadEdge;   // total number of quad edges
   int      numQuadEdgeT;  // total number of quad edges (temp)
   int      numQuadConn;   // number of quadconn elements
   int      numOctFace;    // total number of cells in 3d
   int    **quadConn;      // connectivity data for quadrilaterals
   int     *quad2triList;  // which quad belongs to which original triangle
   int    **quadEdge;      // array for qEdges
   int    **octFace;       // array for oct cell faces

   // array of the grid
   int      numNodePos;    // total number of points
   double  *allNodePos;    // node position of all points in a single array
   double  *surfNodePos;   // node positions of all points on a surface
   double  *cellNodePos;   // node position of all points (if strand)

   // strand grids
   STRAND  *strandGrid;    // strand grids of type strand
   int      numStrandLayer;// total number of strand layers
   int      numConstLayer; // number of constant cell layers (initial)
   double   initMeshLen;   // size of the first cell
   double   meshGrowth;    // mesh spacing for strand grids
   double  *templateDist;  // grid locations on the strand template   

   // mesh quality
   int     *numMeshSkew;   // number of cells of each skewness type
   int    **meshSkewnessID;// skewness of the meshes
   double  *meshSkewness;  // skewness of the meshes
   double  *cellSizeRatio; // ratio of quad area to original triangle”

   // expansion to sub dividing quads
   int      quadLevel;     // division level of quads
   int      nEdgePerSide;  // number of quad edges along a triangular edge
   int      nVertPerSide;  // number of new vertices along a triangular edge
   int      nIntPts;       // additional interior points per triangle
   int     *midEdgeID;     // edge ID for the DO, EO and FO edge for triangle 

   // edge mapping for smoothing at triangular level
   int      numTriQuadEdge;// number of edges at triangular level
   int    **triQuadEdge;   // connectivity for tri edge, mid points and centroid

   // area of triangle and quadrilateral
   double  *triArea;       // area of triangle
   double  *quadArea;      // area of quadrilateral

   // for domain partitioning
   int      boundaryCount; // number of edges sharing a domain
   int     *boundaryID;    // boundary ID for inter-domain boundary

   // for MPI
   int     *ifBoundaryNode; // yes/no for if a node is on a boundary
   int     *boundaryNodeID; // IDs of node on a boundary
   int     *boundaryNodeMap; // grid and node of connecting node on other domain
   int     *revMapNumber;    //
   int      boundaryNodeCount; // number of nodes on subdomain boundary
   double  *boundaryNodePos; // position of nodes on the boundary

   // experimental
   double  *triNodeWeight; // weights for the triangular nodes
   LOOP    *vert1Loop;     // vertex loops for the inner and outer
   LOOP    *vert2Loop;     // loops of a given node (used for smoothing)
   LOOP    *vert3Loop;     // middle loop
   int     *iBoundary;     // index for if boundary node
   int     *q3loop;
   int     *iqloop3;
   int     *flagalter;      

} GRID;

// ==================================================================
#define tracef(x) printf("#meshgen:\t"#x" = %.16e\n",x);
#define trace(x)  printf("#meshgen:\t"#x" = %d\n",x);
#define traces(x) printf("#meshgen:\t"#x"\n");

#define MIN(a,b)  (((a)<(b))?(a):(b));
#define MAX(a,b)  (((a)>(b))?(a):(b));
#define SWAP(a,b) {(a)=(a)+(b);(b)=(a)-(b);(a)=(a)-(b);}
#define SIGN(a)   ((a > 0) - (a < 0));
#endif
// ##################################################################
// END OF FILE
// ##################################################################
