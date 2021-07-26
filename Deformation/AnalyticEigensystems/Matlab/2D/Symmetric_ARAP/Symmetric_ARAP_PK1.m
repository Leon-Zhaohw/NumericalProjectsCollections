function [P] = Symmetric_ARAP_PK1(F)
  [U Sigma V] = svd_rv(F);

  R = U * V';
  S = V * Sigma * V';

  J = det(F);
  J2 = det(F)^2;
  J3 = det(F)^3;
  I1 = trace(S);
  I2 = trace(F' * F);

  P = 2.0 * (1 + 1 / J2) * F - 2 * (1 + 1 / J) * R + (2 / J2) * (I1 - I2 / J) * DJDF(F);
end
