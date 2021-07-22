function showsurface(surface,coordinates,u)
%SHOWSURFACE   Plots surface of three-dimensional body
%    SHOWSURFACE(SURFACE,COORDINATES,U) plots a three-dimensional
%    spline function. SURFACE denotes a set of triangles
%    with dimension (no. of triangles) x 3. This array 
%    include number of nodes. The nodes have to be counted clockwise.
%    or anti-clockwise. Coordinates of nodes are stored in an
%    (no. of coordinates) x 3 - dimensional array called COORDINATES.  
%    Its i'th row defines the x- and y- coordinate. Functionvalues of
%    each node are given by U.

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <showsurface.m> in $(HOME)/acf/fem3d/.

trisurf(surface,coordinates(:,1),coordinates(:,2),coordinates(:,3),u',...
       'facecolor','interp')
axis off
view(160,-30)
