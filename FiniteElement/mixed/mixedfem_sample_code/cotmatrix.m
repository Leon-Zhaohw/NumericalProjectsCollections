function L = cotmatrix(V,F)
% L = cotmatrix(V,F)
% V:  #V x 3 matrix of vertex coordinates
% F:  3 x #F  matrix of indices of triangle corners
% returns  #V x #V matrix of cot weights 
%
% See corresponding paper: "Mixed finite elements for variational surface
% modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
% 2010
%
% Copyright 2010, Alec Jacobson, Denis Zorin, Elif Tosun, NYU
%

  % renaming indices of vertices of triangles for convenience
  i1 = F(1,:); i2 = F(2,:); i3 = F(3,:); 
  % #F x 3 matrices of triangle edge vectors, named after opposite vertices
  v1 = V(i3,:) - V(i2,:);  v2 = V(i1,:) - V(i3,:); v3 = V(i2,:) - V(i1,:);
  % computing areas 
  if size(V,2) == 2
      % 2d vertex data
      dblA = v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1);
  elseif size(V,2) == 3
      %n  = cross(v1,v2,2);  dblA  = multinorm(n,2);

      % area of parallelogram is twice area of triangle
      % area of parallelogram is || v1 x v2 || 
      n  = cross(v1,v2,2); 
      % THIS DOES MATRIX NORM!!! don't use it!!
      % dblA  = norm(n,2);

      % This does correct l2 norm of rows
      dblA = (sqrt(sum((n').^2)))';
  else 
      error('unsupported vertex dimension %d', size(V,2))
  end
  % cotangents and diagonal entries for element matrices
  cot12 = -dot(v1,v2,2)./dblA/2; cot23 = -dot(v2,v3,2)./dblA/2; cot31 = -dot(v3,v1,2)./dblA/2;
  % diag entries computed from the condition that rows of the matrix sum up to 1
  % (follows from  the element matrix formula E_{ij} = (v_i dot v_j)/4/A )
  diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;
  % indices of nonzero elements in the matrix for sparse() constructor
  i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
  j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
  % values corresponding to pairs form (i,j)
  v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
  % for repeated indices (i,j) sparse automatically sums up elements, as we want
  L = sparse(i,j,v,size(V,1),size(V,1));
end
