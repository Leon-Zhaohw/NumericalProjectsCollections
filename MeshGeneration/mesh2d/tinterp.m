function fi = tinterp ( p, t, f, pi, i )

%*****************************************************************************80
%
%% tinterp(): Triangle based linear interpolation.
%
%  Discussion:
%
%    Performs nearest-neighbour extrapolation for points outside the
%    triangulation.
%
%    "dsearch()" replaced by "dsearchn()", John Burkardt, 20 May 2021.
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
%    p  : Nx2 array of nodal XY coordinates, [x1,y1; x2,y2; etc]
%    t  : Mx3 array of triangles as indices, [n11,n12,n13; n21,n22,n23; etc]
%    f  : Nx1 function vector, f(x,y)
%    pi : Jx2 matrix of interpolation points
%
%  Output:
%
%    fi : Jx1 interpolant function vector, fi(xi,yi)
%

%
%  Allocate the output vector.
%
  fi = zeros(size(pi,1),1);
%
%  Deal with points oustide convex hull
%  Do nearest neighbour extrapolation for outside points
%
  out = isnan(i);
  if any(out)
%   d = dsearch(p(:,1),p(:,2),t,pi(out,1),pi(out,2));
    d = dsearchn ( p(:,1:2), t, pi(out,1:2) );
    fi(out) = f(d);
  end
%
%  Keep internal points
%
  pin = pi(~out,:);
  tin = t(i(~out),:);
%
%  Corner nodes
%
  t1 = tin(:,1);
  t2 = tin(:,2);
  t3 = tin(:,3);
%
%  Calculate areas
%
  dp1 = pin-p(t1,:);
  dp2 = pin-p(t2,:);
  dp3 = pin-p(t3,:);
  A3 = abs(dp1(:,1).*dp2(:,2)-dp1(:,2).*dp2(:,1));
  A2 = abs(dp1(:,1).*dp3(:,2)-dp1(:,2).*dp3(:,1));
  A1 = abs(dp3(:,1).*dp2(:,2)-dp3(:,2).*dp2(:,1));
%
%  Linear interpolation
%
  fi(~out) = (A1.*f(t1)+A2.*f(t2)+A3.*f(t3))./(A1+A2+A3);

  return
end      % tinterp()
