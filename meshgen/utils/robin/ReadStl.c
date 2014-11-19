// ##################################################################
//
// ReadStl.m
//
// Code to read STL file and save the data into a file that can
// then be read by a Matlab script
// ##################################################################
#include <stdio.h>
#include <stdlib.h>

int main()
{
   int num_facets;

   // get the number of facets
   getNumFacets(&num_facets);

   // write solution to file that can be read in from Matlab
   writeStlDataToFile(num_facets);

   return 0;
}

// ==================================================================
// Obtain number of facets
// ==================================================================
void getNumFacets(int* numFacets)
{
    // Open surface file
    FILE * fid = fopen("robin_ar2_1160panels.stl","r");
    
    int i = 0;

    fscanf(fid, "%*s %*s");
    do
    {
        fscanf(fid, "%*s %*s %*f %*f %*f"); // facet normal
        fscanf(fid, "%*s %*s");             // outer loop
        fscanf(fid, "%*s %*f %*f %*f");     // vertex 1
        fscanf(fid, "%*s %*f %*f %*f");     // vertex 2
        fscanf(fid, "%*s %*f %*f %*f");     // vertex 3
        fscanf(fid, "%*s");                 // endloop
        fscanf(fid, "%*s");                 // endfacet

        i++;

    } while (fscanf(fid, "%*s %*s") != EOF);

    numFacets[0] = i;

    fclose(fid);
}

// ==================================================================
// Obtain the coordinates and write to file
// ==================================================================
void writeStlDataToFile(int num_facets)
{
    // Open surface file
    FILE * fid1 = fopen("robin_ar2_1160panels.stl","r");
    FILE * fid2 = fopen("dataPoints.dat","w");

    int i = 0;
    double norm_x, norm_y, norm_z;
    double vertex_1_x, vertex_1_y, vertex_1_z;
    double vertex_2_x, vertex_2_y, vertex_2_z;
    double vertex_3_x, vertex_3_y, vertex_3_z;

    // write the total number of data points into file
    fprintf(fid2,"%d %d %d\n",num_facets,0,0);

    // Basic .stl file layout for a facet
    // str str norm(x) norm(y) norm (z)
    // str str
    // str vertex(1,x) vertex(1,y) vertex(1,z)
    // str vertex(2,x) vertex(2,y) vertex(2,z)
    // str vertex(3,x) vertex(3,y) vertex(3,z)
    // str
    // str    
    fscanf(fid1, "%*s %*s");
    do
    {
        fscanf(fid1, "%*s %*s %lf %lf %lf", &norm_x, &norm_y, &norm_z);
        fscanf(fid1, "%*s %*s");
        fscanf(fid1, "%*s %lf %lf %lf", &vertex_1_x, &vertex_1_y, &vertex_1_z);
        fscanf(fid1, "%*s %lf %lf %lf", &vertex_2_x, &vertex_2_y, &vertex_2_z);
        fscanf(fid1, "%*s %lf %lf %lf", &vertex_3_x, &vertex_3_y, &vertex_3_z);
        fscanf(fid1, "%*s");
        fscanf(fid1, "%*s");

        i++;

        // write to file
        fprintf(fid2, "%lf %lf %lf \n", norm_x, norm_y, norm_z);
        fprintf(fid2, "%lf %lf %lf \n", vertex_1_x, vertex_1_y, vertex_1_z);
        fprintf(fid2, "%lf %lf %lf \n", vertex_2_x, vertex_2_y, vertex_2_z);
        fprintf(fid2, "%lf %lf %lf \n", vertex_3_x, vertex_3_y, vertex_3_z);

    } while (fscanf(fid1, "%*s %*s") != EOF);

    fclose(fid1);
    fclose(fid2);

}
// ##################################################################
// END OF FILE
// ##################################################################
