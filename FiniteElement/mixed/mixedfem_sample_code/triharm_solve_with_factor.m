function [xV_solve] = triharm_solve_with_factor( ...
  L, U, P, Q, R, S, M, F, xV_solve, Omega, N0, N1, N2, bndtype, reduction,BZ1,BZ2,V)
  %  Solves Ax = b, where A is factored as PR\AQ = LU and x is the new
  %  coordinates of the vertices, b is determined by S, M and xV using handles
  %  Input:
  %    LUPQR = lu(A)
  %    S = cotmatrix
  %    M = massmatrix
  %    F = faces
  %    xV = x,y,z coordinates of vertices (only exterior is used)
  %    Omega: interior of the domain
  %    N0: boundary of the domain
  %    N1: one layer outside the domain
  %    bndtype:  boundary condition type, ('ext', 'deriv') 
  %      ext means that two rows of x are fixed; deriv means that  1 row of
  %      x is fixed and dx/dn is specified on the same row.
  %    reduction: reduce to a single variable x or keep x,y ('flatten','no_flatten')
  %    BZ1 = bezier vectors (prescribed tangents)
  %    V = original vertex poistions before deformation
  %
  %  Output:
  %    xV = new x,y,z coordinates for all vertices
  %
  %
  % See corresponding paper: "Mixed finite elements for variational surface
  % modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
  % 2010
  %
  % Copyright 2010, Alec Jacobson, NYU
  %

  all = [Omega N0]; 

  % set boundary conditions based on given positions
  % index number is ring number in from interior
  b0 = xV_solve(N0,:);
  f1 = xV_solve(N1,:);
  f2 = xV_solve(N2,:);

  if(strcmp(bndtype,'extxy'))
    % Convert bezier curves into y = Δx = ∂²x/∂t² + ∂²x/∂n²
    % experimental
    size(V)
    size(BZ1)
    size(BZ2)
    y = get_y_from_beziers(F,V,N1,N0,BZ1,BZ2);
  else
    y = zeros(prod(size(N1)),3);
  end

  % solve for each coordinate
  for coord_index = 1:3
    n0 = size(N0,2);
    n1 = size(N1,2);
    n2 = size(N2,2);
    all = [Omega N0]; 
    N12 = [N1 N2];
    n_Omega = size(Omega,2);

    n_all = n_Omega + n0;
    Z_Omega_Omega = sparse(n_Omega, n_Omega);
    Z_all_Omega   = sparse(n_all,   n_Omega);
    Z_all_all     = sparse(n_all,   n_all  );

    if strcmp(bndtype,'deriv')
      error('not implemented');
    else  % bndtype == 'extx'   bndtype == 'extxy' 
      if strcmp(bndtype,'extx')
        y1 = M(N1,N1) \ ( S(N1,N0)*b0(:,coord_index) + ...
          S(N1,N1)*f1(:,coord_index) + S(N1,N2)*f2(:,coord_index));
      else
        y1 = y(:,coord_index);
      end 
      rhs_Dx = -S(all,N0)*b0(:,coord_index) - S(all,N1)*f1(:,coord_index);
      rhs_Dy = -S(all,N1)*y1;
      rhs_Dz = zeros(size(Omega,2),1);%M(Omega,Omega)*g;

      if strcmp(reduction,'flatten')
        rhs = rhs_Dz + S(Omega,all)*( M(all,all) \ ( S(all,all)* ...
          (M(all,all) \ rhs_Dx)  + rhs_Dy));     
      else
          rhs = [ rhs_Dy; rhs_Dz; rhs_Dx];   
      end    
    end    

    % back substitute to solve for this coordinate
    sol = Q*(U\(L\(P*(R\rhs))));

    % extract just primary variable's values 
    if(strcmp(reduction,'no_flatten'))
      if strcmp(bndtype,'deriv') 
        ny = size(Omega,2);
      else
        ny = size(all,2);
      end
      xV_solve(Omega,coord_index) = sol(ny+1:ny+size(Omega,2));   
    else
      xV_solve(Omega,coord_index) = sol(1:end);
    end
  end
end
