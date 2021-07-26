function [P] = Symmetric_Dirichlet_PK1(F)
  [U Sigma V] = svd_rv(F);

  djdf = DJDF(F);
    
  J = det(F);

  C = F' * F;
  IC = trace(C);
  IIC = trace(C * C);
  IIStar = 0.5 * (IC^2 - IIC);
  
  dIIStar = 2 * IC * F - 2 * F * F' * F;
  P = 2 * F + dIIStar / J^2 - (2 / J^3) * IIStar * djdf;
  P = P / 2;
end
