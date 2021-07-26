function [final] = Symmetric_Dirichlet_Psi(F)
  [U Sigma V] = svd_rv(F);
  S = V * Sigma * V';
  R = U * V';

  IC = trace(F' * F);
  J = det(F);

  Finv = F^(-1);
  ICinv = trace(Finv' * Finv);
  final = (1/2) * (IC + ICinv);
end
