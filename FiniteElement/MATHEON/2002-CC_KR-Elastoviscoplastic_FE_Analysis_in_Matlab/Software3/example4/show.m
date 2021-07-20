function show(coordinates,elements,dirichlet,neumann, ...
                Sigma,u,lambda,mu)
AreaOmega = zeros(size(coordinates,1),1);
AvC = zeros(size(coordinates,1),1);
for j = 1:size(elements,1)
  area= det([1,1,1,1;coordinates(elements(j,:),:)'])/6;
  AreaOmega(elements(j,:)) = AreaOmega(elements(j,:)) +area;
  AvC(elements(j,:),:)=AvC(elements(j,:),:)...
  +area*[1;1;1;1]*sqrt(sum(eig(reshape(Sigma(j,:),3,3)).^2));
end;
AvC = AvC./AreaOmega;
E=[dirichlet;neumann]; 
factor=20;
colormap(1-gray);
trisurf(E,factor*u(1:3:size(u,1))+coordinates(:,1), ...
    factor*u(2:3:size(u,1))+coordinates(:,2), ...
    factor*u(3:3:size(u,1))+coordinates(:,3), ...
    AvC, 'facecolor','interp');
 view(-50,35)
% axis equal; axis([-10 70 -100 20 0 40])
 xlabel('x'); ylabel('y'); zlabel('z'); 
%colorbar('vert')




