function M = massmatrix(V,F, type)
  % M = massmatrix(V,F, type)
  % mass matrix for the mesh given by V and F
  % type = 'full': full mass matrix for p.w. linear fem
  % type = 'barycentric': diagonal lumped mass matrix obtained by summing 1/3
  % of areas of surrounding triangles
  % todo: add voronoi-like mass matrix (when the circumcenter is inside
  % triangle, voronoi, otherwise pretend the midpoint of the longest side
  % is the circumcenter)
  %
  % See corresponding paper: "Mixed finite elements for variational surface
  % modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
  % 2010
  %
  % Copyright 2010, Alec Jacobson, Denis Zorin, Elif Tosun, NYU
  %
    
    % renaming indices of vertices of triangles for convenience
    i1 = F(1,:); i2 = F(2,:); i3 = F(3,:); 
    %#F x 3 matrices of triangle edge vectors, named after opposite vertices
    v1 = V(i3,:) - V(i2,:);  v2 = V(i1,:) - V(i3,:); v3 = V(i2,:) - V(i1,:);
    % computing the areas
    if size(V,2) == 2
    % 2d vertex data
      dblA = v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1);
    elseif size(V,2) == 3
      %n  = cross(v1,v2,2);  dblA  = multinorm(n,2);
      n  = cross(v1,v2,2);  
      % THIS DOES MATRIX NORM!!! don't use it!!
      % dblA  = norm(n,2);

      % This does correct l2 norm of rows
      dblA = (sqrt(sum((n').^2)))';
    else 
      error('unsupported vertex dimension %d', size(V,2))
    end
    if strcmp(type,'full')
        % arrays for matrix assembly using 'sparse'
        % indices and values of the element mass matrix entries in the order 
        % (1,2), (2,1),(2,3), (3,2), (3,1), (1,3) (1,1), (2,2), (3,3);
        i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
        j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
        offd_v = dblA/24.;
        diag_v = dblA/12.;
        v = [offd_v,offd_v, offd_v,offd_v, offd_v,offd_v, diag_v,diag_v,diag_v];  
        M = sparse(i,j,v,size(V,1), size(V,1));
    elseif strcmp(type,'barycentric')
        % only diagonal elements
        i = [i1 i2 i3];
        j = [i1 i2 i3];
        diag_v = dblA/6.;
        v = [diag_v,diag_v,diag_v];
        M = sparse(i,j,v,size(V,1), size(V,1));
    elseif strcmp(type,'voronoi')

      % just ported version of intrinsic code

      % edges numbered same as opposite vertices
      FT = F';
      l = [ ...
        sqrt(sum((V(FT(:,2),:)-V(FT(:,3),:)).^2,2)) ...
        sqrt(sum((V(FT(:,3),:)-V(FT(:,1),:)).^2,2)) ...
        sqrt(sum((V(FT(:,1),:)-V(FT(:,2),:)).^2,2)) ...
        ];
      M = massmatrix_intrinsic(l,F,size(V,1),'voronoi');
    else 
        error('bad mass matrix type')
    end
end
