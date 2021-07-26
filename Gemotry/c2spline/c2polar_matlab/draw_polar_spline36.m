% Helper function to display degree (3,6) spline.
function draw_polar_spline36(n, bx, by, bz)

CIRC=[4*n,1:4*n,1:2];
knotsr = [0,0,0,0,1,2,3,4,5,6,7];
knotsc = zeros(1,4*(n+1));
knotsc(1:4:4*(n+1)) = 1:n+1;
knotsc(2:4:4*(n+1)) = 1:n+1;
knotsc(3:4:4*(n+1)) = 1:n+1;
knotsc(4:4:4*(n+1)) = 1:n+1;
knotsc = [0, 0, 0, knotsc, n+2, n+2, n+2];

spx = spmak({knotsr, knotsc}, bx(:,CIRC));
spy = spmak({knotsr, knotsc}, by(:,CIRC));
spz = spmak({knotsr, knotsc}, bz(:,CIRC));

INCR = 0.05;
INCC = 0.02;
xgrid = 0:INCR:4;
ygrid = knotsc(4):INCC:knotsc(end-3);
[X,Y] = meshgrid(xgrid, ygrid);

evx = fnval(spx, [X(:)';Y(:)']);
evy = fnval(spy, [X(:)';Y(:)']);
evz = fnval(spz, [X(:)';Y(:)']);
evx = reshape(evx, size(X));
evy = reshape(evy, size(X));
evz = reshape(evz, size(X));
surf(evx, evy, evz, ones(size(evx)));

shading interp; % Turn on Gouraud shading
hidden off;
lighting gouraud; % Turn on Gouraud shading

view(10, 40);
camlight(-37.5,40);
%camlight(-37.5,-15);
axis off;
title('C2 polar spline - degree (3,6)')

