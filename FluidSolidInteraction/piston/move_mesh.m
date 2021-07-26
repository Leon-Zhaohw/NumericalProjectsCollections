% -----------------------------
% filename move_mesh.m
%
% we compute the mesh deformation

function[wx,vcor_n,vcor_np1]=move_mesh(formulation,Delta_t,du,u_dot_t,vsol,vcor_n,vcor_np1);

[nnt,ndim]=size(vcor_n);
%
% We compute the velocity of each node of the mesh:
%  o Eulerian with mesh deformation technique (formulation = 0)
%  o pure Lagrangian approach (formulation = 1)
%  o ALE for Arbitrary Lagrangian Eulerian (formulation = 2)
%
% For the ALE we use a linear interpolation
%
switch formulation
    case 1
        wx = vsol(:,2)./vsol(:,1);
    case {0,2}
        wx = u_dot_t.*linspace(0,1,nnt)';
end
%
% We update the mesh considering node number 1 as fixed to the wall
% on the left of the chamber.
%
% we first store the coordinates of the previous step
%
vcor_n        = vcor_np1;
%
% next we compute the new coordinates
% applying a linear approximation of the deformation
%
vcor_np1(:,1) = vcor_n(:,1)+wx.*Delta_t;
