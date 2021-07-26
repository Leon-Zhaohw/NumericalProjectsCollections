function [L,U,P,Q,R, S, M] = biharm_factor_system( ...
  V, ...
  F, ...
  bndtype, ...
  masstype, ...
  reduction, ...
  Omega, ...
  N0, ...
  N1)
  % [L,U,P,Q,R, S, M] = biharm_factor_system( ...
  %   V, ...
  %   F, ...
  %   bndtype, ...
  %   masstype, ...
  %   reduction, ...
  %   Omega, ...
  %   N0, ...
  %   N1)
  % construct and factor the system for the biharmonic equation
  % Input,
  % V: #vertices x 3  array of domain positions of vertices (used for stiffness matrix weights) 
  % F: #triangles x 3 array of vertex indices forming triangles
  % bndtype:  boundary condition type, ('ext', 'deriv') 
  %   ext means that two rows of x are fixed; deriv means that  1 row of
  %   x is fixed and dx/dn is specified on the same row.
  % masstype:  type of mass matrix to use ('full', 'barycentric', 'voronoi')
  % reduction: reduce to a single variable x or keep x,y ('flatten','no_flatten')
  % Omega: interior of the domain
  % N0: boundary of the domain
  % N1: one layer outside the domain
  %
  % returns: result of lu factorization of system matrix
  % 
  % notation for indices and vectors from Elif Tosun's thesis
  % 
  % system matrix and rhs will be constructed out of pieces of these
  % matrices   
  %the matrices for the extension and derivative boundary conditions 
  % are slightly different in entries with indices in N0  in the case of 
  % derivative conditions, the entries are computed by integration over the domain Omega only
  % for extension conditions 
  %
  % See corresponding paper: "Mixed finite elements for variational surface
  % modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
  % 2010
  %
  % Copyright 2010, Alec Jacobson, NYU
  %
  interest = [Omega N0 N1 ];
  % Don't both with faces outside of region of interest
  F = F( ...
    ismember(F(:,1),interest) & ...
    ismember(F(:,2),interest) & ...
    ismember(F(:,3),interest),:)';  

  all = [N0 Omega];

  % Build cotangent and mass matrices
  if strcmp(bndtype,'deriv')
    interior_faces = F(:, ismember(F(1,:),all) & ...
                            ismember(F(2,:),all) & ...
                            ismember(F(3,:),all)); 
    S =  cotmatrix(V, interior_faces);  
    M = massmatrix(V, interior_faces, masstype);
  else % bndtype == 'ext'
    S = cotmatrix(V, F);  
    M = massmatrix(V, F, masstype);
  end

  if strcmp(reduction,'flatten')
    % system matrix for the system with x as variable only
    % obtained by eliminating y = M^{-1} (S_{all,Omega} - rhs_Dx)
    A = S(Omega,all) * (M(all,all) \ S(all,Omega)); 
  else % full matrix
    % system matrix for x,y variables
    n_Omega = size(Omega,2);
    Z_Omega_Omega = sparse(n_Omega, n_Omega);
    A = [ -M(all,all)      S(  all,Omega);   ...
           S(Omega,all)    Z_Omega_Omega  ];     
  end

  % factorize A
  [L,U,P,Q,R] = lu(A);
end
