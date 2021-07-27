function [R S] = polar_decomposition_rv(F)
  [U Sigma V] = svd_rv(F);
  R = U * V';
  S = V * Sigma * V';
end
