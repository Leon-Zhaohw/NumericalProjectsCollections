function [xV] = biharm_solve_with_factor( ...
  L, U, P, Q, R, S, M, F, xV, Omega, N0, N1, bndtype, reduction,BZ1,V)
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

  all = [N0 Omega];

  % set boundary conditions based on given positions
  % index number is ring number in from interior
  b0 = xV(N0,:);
  f1 = xV(N1,:);
  if(strcmp(bndtype,'deriv'))
    dxdn = zeros(size(V,1),3);
    bn = [ 2.0*BZ1(N0,:) ; zeros(size(Omega,2),3)];
  else
    bn = zeros(size([N0 Omega],2),3) ;
  end

  % solve for each coordinate
  for coord_index = 1:3
    % build rhs 
    if strcmp(bndtype,'deriv')
      rhs_Dx = -S(all,N0)*b0(:,coord_index) +  bn(:,coord_index); 
    elseif strcmp(bndtype,'ext')
      % in this case, we use the values on the border row
      % to evaluate the rhs
      rhs_Dx = -S(all,N0)*b0(:,coord_index) - S(all,N1)*f1(:,coord_index);
    end
    rhs_Dy = zeros(size(Omega,2),1);
    if strcmp(reduction,'flatten')
        rhs = rhs_Dy + S(Omega,all) * (M(all,all) \ rhs_Dx);
    else % full matrix
        rhs = [ rhs_Dx; rhs_Dy];   
    end
    % back substitute to solve for this coordinate
    sol = Q*(U\(L\(P*(R\rhs))));

    % extract just primary variable's values 
    if(strcmp(reduction,'no_flatten'))
      ny = size(all,2);
      xV(Omega,coord_index) = sol(ny+1:end);   
    else
      xV(Omega,coord_index) = sol(1:end);
    end
  end
end
