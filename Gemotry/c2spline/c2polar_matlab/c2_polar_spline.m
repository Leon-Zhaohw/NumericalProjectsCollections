% Note, this function only works on a single coordinate of the polar
% configuration (i.e. x, y, or z).
%
% Input: 2-layer polar configuration polcfg with 4 x n vertices,
%        where polcfg(1,:) is the same vertex
% Output: Spline control net b of C^2 degree (3, 6) spline with knots
%        [0,0,0,0,1,2,3,4,5,6,7] and
%        [0,0,0,0, 1,1,1,1, 2,2,2,2, ..., n-1,n-1,n-1,n-1] (periodic)
%        and 7 x 4n control points.
%        b(1,:) is the collapsed boundary that corresponds to the pole.
function b = c2_polar_spline(n, polcfg)

c = cos(2*pi/n * [0:n-1]);
s = sin(2*pi/n * [0:n-1]);
c2 = cos(4*pi/n * [0:n-1]);
s2 = sin(4*pi/n * [0:n-1]);
one = ones(1,n);

q = polar_subdiv(n, polcfg);

p0 = (2.0/3.0)*q(1,1) + (1.0/3.0)/n*sum(q(2,:));
p1 = (2.0/n) * c * q(2,:)';
p2 = (2.0/n) * s * q(2,:)';
p3 = -q(1,1) + 1/n*sum(q(2,:));
p4 = (2.0/n) * c2 * q(2,:)';
p5 = (2.0/n) * s2 * q(2,:)';

v0 = bsp_mul33(one, one);
v1 = bsp_mul33(c, one);
v2 = bsp_mul33(s, one);
v3 = bsp_mul33(c, c) + bsp_mul33(s, s);
v4 = bsp_mul33(c, c) - bsp_mul33(s, s);
v5 = 2*bsp_mul33(c, s);

b = zeros(7,4*n);
b(1,:) = p0*v0;
b(2,:) = p0*v0 + (1/3)*(p1*v1 + p2*v2);
b(3,:) = p0*v0 +       (p1*v1 + p2*v2) + (2/3)*(p3*v3 + p4*v4 + p5*v5);
b(4,:) = bsp_mul33(one, q(3,:));
b(5,:) = bsp_mul33(one, q(4,:));
b(6,:) = bsp_mul33(one, q(5,:));
b(7,:) = bsp_mul33(one, q(6,:));

