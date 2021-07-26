function [final] = ARAP_Psi(F)
  [U Sigma V] = svd_rv(F);
  S = V * Sigma * V';
  R = U * V';
  final = norm(F - R, 'fro')^2;
end