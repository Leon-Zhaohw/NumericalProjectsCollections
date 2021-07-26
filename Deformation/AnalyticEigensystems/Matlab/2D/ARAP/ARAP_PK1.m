function [P] = ARAP_PK1(F)
  [U Sigma V] = svd_rv(F);
  S = V * Sigma * V';
  R = U * V';
  
  P = R * (2 * (S - eye(2,2)));
end
