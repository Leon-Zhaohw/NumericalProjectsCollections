The included code is for a 2D implementation of "Matching Fluid Simulation Elements to Surface Geometry and Topology" from SIGGRAPH 2010.

Various simulation parameters can be tweaked in "main.cpp" (including a selection of input geometries as signed distance fields, laid out in "scenes.cpp".)

The code is released into the public domain, with the exception of Jonathan Shewchuk's Triangle 2D meshing library (which is used for Delaunay triangulation). See triangle.c for Triangle's licensing information.

WINDOWS
---------

Project files for Visual Studio 2008 are included, however you will need to modify it to link to a BLAS library (or write your own implementations of a few required functions). Any code you might need to change in that regard is located in "blaswrapper.h".

OS/X AND LINUX
---------

To build on OS/X or Linux using make, you will need to create a Makefile.local_defs file.  See the included Makefile.example_defs as a guide.  Then, to build the binary executable, build the dependencies and the main program.  E.g.:

$ make depend
$ make release

This should produce an optimized executable called voronoi2d_release.  You can also build an unoptimized executable with debug symobls by running "make debug".

---------

Please report any bugs or improvements to "tbrochu@cs.ubc.ca" or "christopherbatty@yahoo.com".


