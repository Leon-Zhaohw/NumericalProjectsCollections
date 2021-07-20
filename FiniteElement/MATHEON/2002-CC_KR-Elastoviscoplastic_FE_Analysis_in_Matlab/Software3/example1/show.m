function show(coordinates,elements,Sigma,u,lambda,mu)
AreaOmega = zeros(size(coordinates,1),1);
AvS = zeros(size(coordinates,1),4);
for j = 1:size(elements,1)
  area= det([1,1,1;coordinates(elements(j,:),:)'])/2;
  AreaOmega(elements(j,:)) = AreaOmega(elements(j,:)) +area;
  AvS(elements(j,:),:) = AvS(elements(j,:),:) +area*[1;1;1]*Sigma(j,:);
end;
AvS = AvS./(AreaOmega*[1,1,1,1]);
AvC=(mu/(24*(mu+lambda)^2)+1/(8*mu))*(AvS(:,1)+...
     AvS(:,4)).^2+1/(2*mu)*(AvS(:,2).^2-AvS(:,1).*AvS(:,4));
factor=20;
colormap(1-gray);
trisurf(elements,factor*u(1:2:size(u,1))+coordinates(:,1), ...
    factor*u(2:2:size(u,1))+coordinates(:,2), ...
    zeros(size(coordinates,1),1), AvC, 'facecolor','interp');
view(0,90);
colorbar('vert')




