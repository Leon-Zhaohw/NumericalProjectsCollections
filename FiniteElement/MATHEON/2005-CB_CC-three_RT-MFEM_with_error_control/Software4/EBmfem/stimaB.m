function B=stimaB(coord);
%N=coord(:)*ones(1,3)-[coord;coord;coord];
N=coord(:)*ones(1,3)-repmat(coord,3,1);
D=diag([norm(N([5,6],2)) norm(N([1,2],3)) norm(N([1,2],2))]);
M=spdiags([ones(6,1),ones(6,1),2*ones(6,1),ones(6,1),ones(6,1)],...
          [-4,-2,0,2,4],6,6);
B = D*N'*M*N*D/(24*det([1,1,1;coord])); 