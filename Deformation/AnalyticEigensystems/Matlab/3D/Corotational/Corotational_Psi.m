function [final] = Corotational_Psi(F)
  [U Sigma V] = svd_rv(F);
  S = V * Sigma * V';
  R = U * V';
  final = norm(F - R, 'fro')^2 + (1/2) * trace(S - eye(3,3))^2;
end