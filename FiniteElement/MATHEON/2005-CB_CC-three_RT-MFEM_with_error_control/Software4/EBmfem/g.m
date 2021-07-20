function val = g(x,n);
a = angle((x(:,1)+x(:,2)*i)*(-1-i)/sqrt(2))+pi*3/4;
r = sqrt(x(:,1).^2+x(:,2).^2);
val=(2/3*r.^(-1/3).*[-sin(a/3),cos(a/3)])*n';
  