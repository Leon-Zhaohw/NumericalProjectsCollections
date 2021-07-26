% Helper function to display degree (3,3) spline.
function draw_polar_spline33(n, bx, by, bz)

CIRC=[n,1:n,1:2];
knotsr = [0,0,0,0,1,2,3,4,5];
knotsc = 1:length(CIRC)+4;

spx = spmak({knotsr, knotsc}, bx(:,CIRC));
spy = spmak({knotsr, knotsc}, by(:,CIRC));
spz = spmak({knotsr, knotsc}, bz(:,CIRC));

INCR = 0.05;
INCC = 0.02;
xgrid = 0:INCR:2;
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
title('C1 polar spline - degree (3,3)')

