
function u_h=tangent(coordinate,dirichlet,u_h);

Tangent = zeros(2*size(coordinate,1),2);
eps_in_direction_edge_d = zeros(2*size(coordinate,1),2);
 
for j = 1 : size(dirichlet,1)
   TangentToEdge = (coordinate(dirichlet(j,2),:)-coordinate(dirichlet(j,1),:)) * ... %CB   
             norm(coordinate(dirichlet(j,2),:)-coordinate(dirichlet(j,1),:));
   tan_eps =  (1/1000000) * (coordinate(dirichlet(j,2),:)-coordinate(dirichlet(j,1),:)) /...
            norm(coordinate(dirichlet(j,2),:)-coordinate(dirichlet(j,1),:));
   if norm(Tangent(2*dirichlet(j,1)-1,:)) > 0
     Tangent(2*dirichlet(j,1),:) = TangentToEdge;
      eps_in_direction_edge_d(2*dirichlet(j,1),:) = tan_eps;
  else
     Tangent(2*dirichlet(j,1)-1,:) = TangentToEdge;
     eps_in_direction_edge_d(2*dirichlet(j,1)-1,:) = tan_eps;
 end
   if norm(Tangent(2*dirichlet(j,2)-1,:)) > 0
     Tangent(2*dirichlet(j,2),:) = TangentToEdge;
     eps_in_direction_edge_d(2*dirichlet(j,2),:) = - tan_eps;
 else
     Tangent(2*dirichlet(j,2)-1,:) = TangentToEdge;
     eps_in_direction_edge_d(2*dirichlet(j,2)-1,:) = - tan_eps; 
 end
end
 maske = zeros(size(coordinate,1),1);
 maske(dirichlet) = ones(size(dirichlet));
 dirichletNode = find(maske);
 for j = 1 : size(dirichletNode,1)
   if det(Tangent([2*dirichletNode(j)-1:2*dirichletNode(j)],:)) ~= 0
       u_h(dirichletNode(j),:) = (Tangent([2*dirichletNode(j)-1:2*dirichletNode(j)],:) \ ...
           [duds(coordinate(dirichletNode(j),:) + eps_in_direction_edge_d(2*dirichletNode(j)-1,:));...
            duds(coordinate(dirichletNode(j),:)+eps_in_direction_edge_d(2*dirichletNode(j),:)) ])';
   else 
       normal = Tangent(2*dirichletNode(j)-1,:) * [0 1;-1 0];
       u_h(dirichletNode(j),:) = ([Tangent(2*dirichletNode(j)-1,:);normal] \ ...
          [duds(coordinate(dirichletNode(j),:));u_h(dirichletNode(j),:) * normal'])'; 
   end
 end
 