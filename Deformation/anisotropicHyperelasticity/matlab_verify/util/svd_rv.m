% Rotation-variant SVD which pulls any reflections out of
% U and V and pushes them to Sigma
function [U Sigma V] = svd_rv(F)
  [U Sigma V] = svd(F);
  
  % reflection matrix
  L = eye(3,3);
  L(3,3) = det(U * V');
  
  % see where to pull the reflection out of
  detU = det(U);
  detV = det(V);
  if (detU < 0 && detV > 0)
      U = U * L;
  elseif (detU > 0 && detV < 0)
      V = V * L;
  elseif (detU < 0 && detV < 0)
      L = [1 0 0; 
           0 1 0;
           0 0 -1];
      U = U * L;
      V = V * L;
      L = eye(3,3);
  end  
  
  % push the reflection to the diagonal
  Sigma = Sigma * L;
end
