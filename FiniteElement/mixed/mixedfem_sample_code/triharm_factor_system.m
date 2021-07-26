function [L,U,P,Q,R, S, M] = triharm_factor_system( ...
  V, ...
  F, ...
  bndtype, ...
  masstype, ...
  reduction, ...
  Omega, ...
  N0, ...
  N1, ...
  N2)
  % [L,U,P,Q,R, S, M] = triharm_factor_system( ...
  %   V, ...
  %   F, ...
  %   bndtype, ...
  %   masstype, ...
  %   reduction, ...
  %   Omega, ...
  %   N0, ...
  %   N1, ...
  %   N2)
  % construct and factor the system for the biharmonic equation
  % Input,
  % V: #vertices x 3  array of domain positions of vertices (used for stiffness matrix weights) 
  % F: #triangles x 3 array of vertex indices forming triangles
  % bndtype:  boundary condition type, ('ext', 'deriv') 
  %   extxy means that three rows of x are fixed; extxy means that 2 rows are
  %   fixed
  %   x is fixed and d2x/dn2 is specified on the same row.
  % masstype:  type of mass matrix to use ('full', 'barycentric', 'voronoi')
  % reduction: reduce to a single variable x or keep x,y ('flatten','no_flatten')
  % Omega: interior of the domain
  % N0: boundary of the domain
  % N1: one layer outside the domain
  % N2: two layers outside the domain
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

  interest = [Omega N0 N1 N2];
  % Don't both with faces outside of region of interest
  F = F( ...
    ismember(F(:,1),interest) & ...
    ismember(F(:,2),interest) & ...
    ismember(F(:,3),interest),:)';  

  n0 = size(N0,2);
  n1 = size(N1,2);
  n2 = size(N2,2);
  n_Omega = size(Omega,2);
  all = [Omega N0];
  N12 = [N1 N2];

  % Build cotangent and mass matrices
  if strcmp(bndtype,'deriv')
    fprintf('Probably isnt working...\n');
    interior_faces = F(:, ismember(F(1,:),all) & ...
                            ismember(F(2,:),all) & ...
                            ismember(F(3,:),all)); 
    S =  cotmatrix(V, interior_faces);  
    M = massmatrix(V, interior_faces, masstype);
  else % bndtype == 'extx' || bndtype == 'extxy'
    S = cotmatrix(V, F);  
    M = massmatrix(V, F, masstype);
  end

  n_all = n_Omega + n0;
  Z_Omega_Omega = sparse(n_Omega, n_Omega);
  Z_all_Omega   = sparse(n_all,   n_Omega);
  Z_all_all     = sparse(n_all,   n_all  );

  if  strcmp(bndtype,'deriv')
    if strcmp(reduction,'flatten')
      error('not implemented');
    else % no flatten
     A = [S(Omega,Omega) Z_Omega_Omega   -M(Omega,all);    ...
       Z_Omega_Omega  Z_Omega_Omega    S(Omega,all);   ...
       -M(all,Omega)    S(all,Omega)   Z_all_all   ];     
    end
  else  % bndtype == 'extx'   bndtype == 'extxy' 
    if strcmp(reduction,'flatten')
      % system matrix for the system with x as variable only
      % obtained by eliminating y,z 
      A = S(Omega,all) * (M(all,all) \ ( S(all,all)* (M(all,all) \ S(all, Omega))));
    else % full matrix
       % system matrix for x,y variables
      A = [S(all,all)    Z_all_Omega   -M(all,all);    ...
        Z_all_Omega'  Z_Omega_Omega  S(Omega,all);   ...
        -M(all,all)    S(all,Omega)   Z_all_all   ];
    end    
  end    

  % factorize A
  [L,U,P,Q,R] = lu(A);
end
