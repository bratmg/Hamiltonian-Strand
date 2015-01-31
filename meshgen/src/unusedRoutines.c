// ################################################################## 
//
// Contains a bunch of routines that are not used
//
//
// ################################################################## 

// ################################################################## 
// computeNodeWeights (originally from smoothOperations.c)
// 
// Computes the weights associated with each node
// ##################################################################
void computeNodeWeights(GRID *g)
{
   printf("#meshgen: Computing node weights ...\n");
   // logic - each node - find triangle loops - 
   //         each triangle - find edges
   //         if node part of an edge - compute associated length
   //         rinse and repeat
   int    i,j,k,i1,i2,nv;
   int    vertID;
   double posA[3],posB[3],sum;

   g->triNodeWeight = (double *) malloc(sizeof(double)*g->numTriNode);
   
   // loop through all nodes
   for (i = 0; i < g->numTriNode; i++)
   {
      posA[0] = g->nodePosTri[3*i  ];
      posA[1] = g->nodePosTri[3*i+1];
      posA[2] = g->nodePosTri[3*i+2];

      sum = 0.0;
      i1 = g->vertLoop->index[i];
      i2 = g->vertLoop->index[i+1];
      nv = i2-i1;
      // loop through the vertices of a particular node
      for (j = i1; j < i2; j++)
      {
         // ID of the vertex
         vertID = g->vertLoop->ID[j];

         // position of vertex
         posB[0] = g->nodePosTri[3*vertID  ];
         posB[1] = g->nodePosTri[3*vertID+1];
         posB[2] = g->nodePosTri[3*vertID+2];

         sum += (posA[0]-posB[0])*(posA[0]-posB[0]) +
                (posA[1]-posB[1])*(posA[1]-posB[1]) +
                (posA[2]-posB[2])*(posA[2]-posB[2]);
      }
      sum = sqrt(sum)/nv;
      g->triNodeWeight[i] = sum;

   } // i loop


}


// ##################################################################
// Originally from smoothLoops.c
//
// smooth the loops and not the nodes
//
// ##################################################################
void smoothLoops(GRID *g)
{

   printf("#meshgen: Smoothing along loops ...\n");

   int     i,ii,j,k1,k2,m,kk;
   int     msweep,temp0,temp1,temp2,threeNum;
   int     i1,i2,nel,e1,e2,is,ie;
   int     edge1[2],edge2[2];
   int     *vert4ID,*vert4Index;
   double *x1,sum1,sum2,sum3;
   x1 = (double *) malloc(sizeof(double)*3*g->numNodePos);
   
   threeNum = 3*g->numNodePos;
   for (i = 0; i < threeNum; i++)
      x1[i] = g->allNodePos[i];

   g->vert1Loop = (LOOP *) malloc(sizeof(LOOP));
   g->vert2Loop = (LOOP *) malloc(sizeof(LOOP));
   g->vert3Loop = (LOOP *) malloc(sizeof(LOOP));

   // set the total length of these loops
   g->vert1Loop->totLen = 0.5*(g->quadLoop->totLen);
   g->vert2Loop->totLen = 0.5*(g->quadLoop->totLen);
   g->vert3Loop->totLen = 0.5*(g->quadLoop->totLen);

   g->vert1Loop->ID = (int *) malloc(sizeof(int)*g->vert1Loop->totLen);
   g->vert2Loop->ID = (int *) malloc(sizeof(int)*g->vert2Loop->totLen);
   g->vert3Loop->ID = (int *) malloc(sizeof(int)*g->vert3Loop->totLen);
   vert4ID          = (int *) malloc(sizeof(int)*g->vert3Loop->totLen);

   g->vert1Loop->index = (int *) malloc(sizeof(int)*(1+g->numTriNode));
   g->vert2Loop->index = (int *) malloc(sizeof(int)*(1+g->numTriNode));
   g->vert3Loop->index = (int *) malloc(sizeof(int)*(1+g->numTriNode));
   vert4Index          = (int *) malloc(sizeof(int)*(1+g->numTriNode));

   int iouter  = 1;
   int iinner  = 1;
   int imiddle = 1;
   int iinnest = 1;

   // Initialize
   for (i = 0; i <= g->numTriNode; i++)
   {
      g->vert1Loop->index[i] = g->vert2Loop->index[i] = 0;
      g->vert3Loop->index[i] = vert4Index[i] = 0;
   }

// ==================================================================
// Set up the indices for the different loops
// 4 loops per node
// ==================================================================   
   k1 = 0; k2 = 0; kk = 0;
   // Loop through the number of nodes
   for (i = 0; i < g->numTriNode; i++)
   {
      //ii = g->colourIndex[i];

      // ============================================================
      // INNERMOST LOOP
      // ============================================================
      // indices for start and end of the inner and outer loop
      // for a given triangular node
      i1  = g->iqloop3[i];
      i2  = g->iqloop3[i+1];      
      nel = i2-i1;

      // Loop through the inner
      for (j = i1; j < i2; j++)
      {
         e1 = g->q3loop[j];         

         edge1[0] = g->quadEdge[e1][0];
         edge1[1] = g->quadEdge[e1][1];

         //vert4ID[k2] = edge1[0];
         vert4ID[k2] = edge1[1];
         k2++;
      } // j loop
      vert4Index[kk+1] = vert4Index[kk] + nel;
   
      // ============================================================
      // 3 other loops
      // ============================================================
      // indices for start and end of the inner and outer loop
      // for a given triangular node
      i1  = g->quadLoop->index[2*i  ];
      i2  = g->quadLoop->index[2*i+1];            
      nel = i2-i1;      

      // Loop through the inner
      for (j = i1; j < i2; j++)
      {
         e1 = g->quadLoop->ID[j];
         e2 = g->quadLoop->ID[j+nel];

         edge1[0] = g->quadEdge[e1][0];
         edge1[1] = g->quadEdge[e1][1];

         edge2[0] = g->quadEdge[e2][0];
         edge2[1] = g->quadEdge[e2][1];

         // Note that edge1[1] and edge2[0] are always
         // the same. Therefore the other remaining
         // vertex IDs are those of the two loops
         g->vert1Loop->ID[k1] = edge1[0];
         g->vert2Loop->ID[k1] = edge2[1];
         g->vert3Loop->ID[k1] = edge1[1];
         k1++;
      } // j loop


      g->vert1Loop->index[kk+1] = g->vert1Loop->index[kk] + nel;
      g->vert2Loop->index[kk+1] = g->vert2Loop->index[kk] + nel;
      g->vert3Loop->index[kk+1] = g->vert3Loop->index[kk] + nel;

      
      kk++;


   } // i loop

// ==================================================================
// With the created vert loops - apply laplacian smoothing
// to each of the loops
// Perform the sweeps   
// ==================================================================   
   msweep = 5;
   traces('are you sure that msweep is 5?');
   exit(1);
   for (m = 0; m < msweep; m++)
   {
      for (i = 0; i < g->numTriNode; i++)
      {         
         
         // =========================================================
         // INNERMOST LOOP
         // =========================================================
         i1 = vert4Index[i];
         i2 = vert4Index[i+1];
         nel= i2-i1;
         sum1=sum2=sum3=0.0;nel=0;
         // smooth the positions (same for open or closed loop)
         for (j = i1+1; j < i2-2; j++)
         {  
            temp0 = 3*vert4ID[j];
            temp1 = 3*vert4ID[j-1];
            temp2 = 3*vert4ID[j+1];
            
            if(iinnest==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
               sum1 += g->allNodePos[temp0];
               sum2 += g->allNodePos[temp0+1];
               sum3 += g->allNodePos[temp0+2];
               nel++;
            }

         } // j loop
         
         if( vert4ID[i1] == vert4ID[i2-1])
         {
            temp0 = 3*vert4ID[i1];
            temp1 = 3*vert4ID[i1+1];
            temp2 = 3*vert4ID[i2-1];

            if(iinnest==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
               sum1 += g->allNodePos[temp0];
               sum2 += g->allNodePos[temp0+1];
               sum3 += g->allNodePos[temp0+2];
               nel++;
            }
         
         }
         ii = g->colourIndex[i];
         if(iinnest==1)
         {
            x1[3*ii  ] = sum1/nel;
            x1[3*ii+1] = sum2/nel;
            x1[3*ii+2] = sum3/nel;
         }
         
         // =========================================================
         // OUTER 3 LOOPS
         // =========================================================
         i1  = g->vert1Loop->index[i];
         i2  = g->vert1Loop->index[i+1];
         nel = i2-i1;

         // smooth the positions (same for open or closed loop)
         for (j = i1+1; j < i2-2; j++)
         {  
            
            // OUTER LOOP       
            temp0 = 3*g->vert1Loop->ID[j];
            temp1 = 3*g->vert1Loop->ID[j-1];
            temp2 = 3*g->vert1Loop->ID[j+1];
            
            if(iouter==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }
            // INNER LOOP         
            temp0 = 3*g->vert2Loop->ID[j];
            temp1 = 3*g->vert2Loop->ID[j-1];
            temp2 = 3*g->vert2Loop->ID[j+1];
            
            if(iinner==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }

            // MIDDLE LOOP
            temp0 = 3*g->vert3Loop->ID[j];
            temp1 = 3*g->vert3Loop->ID[j-1];
            temp2 = 3*g->vert3Loop->ID[j+1];
            
            if(imiddle==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }
         } // j loop
         
         // check if closed loop or not. Done by checking the
         // vert ID of the beginning and end of the loop. Can use
         // wither vert1Loop or vert2Loop      
         
         if( g->vert1Loop->ID[i1] == g->vert1Loop->ID[i2-1])
         {
            temp0 = 3*g->vert1Loop->ID[i1];
            temp1 = 3*g->vert1Loop->ID[i1+1];
            temp2 = 3*g->vert1Loop->ID[i2-1];

            if(iouter==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }

            // INNER LOOP
            temp0 = 3*g->vert2Loop->ID[i1];
            temp1 = 3*g->vert2Loop->ID[i1+1];
            temp2 = 3*g->vert2Loop->ID[i2-1];
            
            if(iinner==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }

            // MIDDLE LOOP
            temp0 = 3*g->vert3Loop->ID[i1];
            temp1 = 3*g->vert3Loop->ID[i1+1];
            temp2 = 3*g->vert3Loop->ID[i2-1];
            
            if(imiddle==1)
            {
               x1[temp0  ] = 0.5*(g->allNodePos[temp1  ] + g->allNodePos[temp2  ]);
               x1[temp0+1] = 0.5*(g->allNodePos[temp1+1] + g->allNodePos[temp2+1]);
               x1[temp0+2] = 0.5*(g->allNodePos[temp1+2] + g->allNodePos[temp2+2]);
            }
         }

      } // i loop

// ==================================================================         
// POSITION UPDATE
// ==================================================================            
      for (i = 0; i < g->numNodePos; i++)
      {         
         //i = g->colourIndex[ii];
         g->allNodePos[3*i  ] = g->iBoundary[i]*g->allNodePos[3*i  ]+x1[3*i  ]*(1-g->iBoundary[i]);
         g->allNodePos[3*i+1] = g->iBoundary[i]*g->allNodePos[3*i+1]+x1[3*i+1]*(1-g->iBoundary[i]);
         g->allNodePos[3*i+2] = g->iBoundary[i]*g->allNodePos[3*i+2]+x1[3*i+2]*(1-g->iBoundary[i]);
      }
      
   } // m loop

}

