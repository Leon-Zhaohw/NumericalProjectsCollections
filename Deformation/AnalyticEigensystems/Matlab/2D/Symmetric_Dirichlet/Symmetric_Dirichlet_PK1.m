function [P] = Symmetric_Dirichlet_PK1(F)
  [U Sigma V] = svd_rv(F);
  J = det(F);
  J3 = det(F)^3;
  I2 = trace(F' * F);

  P = (1 + 1 / J^2) * F - I2 / J3 * DJDF(F);
end
