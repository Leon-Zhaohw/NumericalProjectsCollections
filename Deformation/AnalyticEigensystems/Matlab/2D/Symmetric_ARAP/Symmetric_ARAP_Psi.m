function [final] = Symmetric_ARAP_Psi(F)
  [U Sigma V] = svd_rv(F);
  S = V * Sigma * V';
  R = U * V';

  C = F' * F;
  IC = trace(F' * F);
  IIC = trace(C' * C);
  IIStarC = 0.5 * (IC^2 - IIC);
  J = det(F);
  IS = trace(S);
  djdf = DJDF(F);

  FmR = F - R;
  FinvmR = F^-1 - R';
  final = norm(FmR, 'fro')^2 + norm(FinvmR, 'fro')^2;
end
