function A = triarea(p,t)

%*****************************************************************************80
%
%% triarea(): Area of triangles assuming counter-clockwise (CCW) node ordering.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    13 March 2013
%
%  Author:
%
%    Darren Engwirda
%
%  Input:
%
%  P  : Nx2 array of XY node co-ordinates
%  T  : Mx3 array of triangles as indices into P
%
%  Output:
%
%  A  : Mx1 array of triangle areas
%
  d12 = p(t(:,2),:)-p(t(:,1),:);
  d13 = p(t(:,3),:)-p(t(:,1),:);
  A = 0.5 * (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1));

  return
end
