% Note, this function only works on a single coordinate of the polar
% configuration (i.e. x, y, or z).
%
% Input: 2-layer polar configuration polcfg with 4 x n vertices,
%        where polcfg(1,:) is the same vertex
% Output: Spline control net b of C^1 bicubic spline with knots
%        [0,0,0,0,1,2,3,4] and
%        [0,1,2, ..., n-1] (periodic)
%        and 5 x n control points.
%        b(1,:) is the collapsed boundary that corresponds to the pole.
function b = c1_polar_spline(n, polcfg)

q = polcfg;
c = cos(2*pi/n * [0:n-1]);
s = sin(2*pi/n * [0:n-1]);
p0 = (2.0/3.0)*q(1,1) + (1.0/3.0)/n*sum(q(2,:));
p1 = (2.0/n) * c * q(2,:)';
p2 = (2.0/n) * s * q(2,:)';

v0 = ones(1,n);
v1 = c;
v2 = s;

b = zeros(5,n);
b(1,:) = p0 * v0;
b(2,:) = p0 * v0 + (1/3)*(p1*v1 + p2*v2);
b([3,4,5],:) = q([2,3,4],:);

