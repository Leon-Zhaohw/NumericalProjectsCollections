This is an implementation of the technique proposed in
    Myles, Ni, Peters. "Fast Parallel Construction of Smooth Surfaces from
    Meshes with Tri/Quad/Pent Facets", SGP 2008.

This code has been tested to work properly with
- Windows Vista 32- and 64-bits
- Visual Studio 2005
- November 2007 DirectX SDK
- CGAL 3.3.1
- GeForce 8800GTX, 8800GT (used for the paper), and ATI Radeon HD 4870 X2

NOTES:
- The performance information on the original paper used Gouraud (smooth)
  shading, whereas Phong shading is on by default in this code. This can be
  changed by modifying the line "#define PHONG 1" at the beginning of
  surface.fx.
