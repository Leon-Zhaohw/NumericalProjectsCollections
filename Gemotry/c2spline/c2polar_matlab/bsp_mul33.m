% Helper function for c2_polar_spline.m

% Multiplies two uniform periodic cubic splines each with n control points
% to get a periodic spline of degree 6 with 4-fold knots.
function c = bsp_mul33(a, b)

n = length(a);

a0 = a([end,1:end-1]);
a1 = a;
a2 = a([2:end,1]);
a3 = a([3:end,1,2]);

b0 = b([end,1:end-1]);
b1 = b;
b2 = b([2:end,1]);
b3 = b([3:end,1,2]);

c = zeros(1,4*n);
c(1:4:end) = a0 .* b1 / 10.0 + a0 .* b2 / 30.0 + a1 .* b0 / 10.0 + 8.0 / 15.0 * a1 .* b1 + a1 .* b2 / 10.0 + a2 .* b0 / 30.0 + a2 .* b1 / 10.0;
c(2:4:end) = a0 .* b1 / 90.0 + a0 .* b2 / 45.0 + 16.0 / 45.0 * a1 .* b1 + 7.0 / 30.0 * a1 .* b2 + 7.0 / 30.0 * a2 .* b1 + a2 .* b2 / 9.0 + a1 .* b0 / 90.0 + a2 .* b0 / 45.0;
c(3:4:end) = a0 .* b3 / 720.0 + a1 .* b3 / 180.0 + a2 .* b3 / 720.0 + a0 .* b1 / 720.0 + a0 .* b2 / 180.0 + a1 .* b0 / 720.0 + 19.0 / 90.0 * a1 .* b1 + 197.0 / 720.0 * a1 .* b2 + a2 .* b0 / 180.0 + 197.0 / 720.0 * a2 .* b1 + 19.0 / 90.0 * a2 .* b2 + a3 .* b0 / 720.0 + a3 .* b1 / 180.0 + a3 .* b2 / 720.0;
c(4:4:end) = a1 .* b1 / 9.0 + 7.0 / 30.0 * a1 .* b2 + a1 .* b3 / 45.0 + 7.0 / 30.0 * a2 .* b1 + 16.0 / 45.0 * a2 .* b2 + a2 .* b3 / 90.0 + a3 .* b1 / 45.0 + a3 .* b2 / 90.0;

