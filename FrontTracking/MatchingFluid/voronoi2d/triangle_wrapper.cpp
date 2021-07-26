#include "triangle_wrapper.h"

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <stdio.h>
#include <stdlib.h>

extern "C" {
   #include "triangle/triangle.h"
}

void compute_delaunay(const std::vector<Vec2f>& points, std::vector<Vec3ui>& tris) {

   triangulateio in, out;

   in.numberofpoints = points.size();
   in.numberofpointattributes = 0;
   in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
   for(unsigned int i = 0; i < points.size(); ++i) {
      in.pointlist[2*i] = points[i][0];
      in.pointlist[2*i+1] = points[i][1];
   }

   in.pointmarkerlist = (int *) NULL;

   in.numberofsegments = 0;
   in.numberofholes = 0;
   in.numberofregions = 0;

   out.trianglelist = (int *) NULL;
   out.triangleattributelist = (REAL *) NULL;

   //Delaunay triangle the point set, without generating the vertices (N)
   //or boundary markers (B), dumping lots of text output(Q),
   //or inserting Steiner points (Y).
   //Also, number things starting from zero (z).

   char options[256];
#ifdef WIN32
   _snprintf( options, 256, "YBNQz" );
#else
   snprintf( options, 256, "YBNQz" );
#endif
   
   triangulate( options, &in, &out, (triangulateio*)NULL);

   for(int i = 0; i < out.numberoftriangles; ++i) {      
      Vec3ui tri(out.trianglelist[i*3], out.trianglelist[i*3+1], out.trianglelist[i*3+2]);
      tris.push_back(tri);
   }
   
   //free everything, including stuff allocated by Triangle.
   //I think these are the only two to worry about.
   free(in.pointlist);
   free(out.trianglelist);
   

}