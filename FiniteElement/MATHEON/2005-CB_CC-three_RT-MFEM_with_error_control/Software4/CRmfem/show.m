function ShowDisplacement(element,coordinate,nodes2edge,x)
 X = reshape(coordinate(element',1), size(element,2), size(element,1));
 Y = reshape(coordinate(element',2), size(element,2), size(element,1));
 Z = zeros(size(element,2), size(element,1));
 for j = 1:size(element,1)
   I = diag(nodes2edge(element(j,[2 3 1]), element(j,[3 1 2])));
   Z(:,j) = [-1, 1, 1; 1, -1, 1; 1, 1, -1]*x(I);
 end;
 fill3(X,Y,Z,'w');
 title('Solution u')
 view(150,40);