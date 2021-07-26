function [P] = Symmetric_ARAP_PK1(F)
  [U Sigma V] = svd_rv(F);
  R = U * V';
  S = V * Sigma * V';

  djdf = DJDF(F);
  J = det(F);

  C = F' * F;
  IC = trace(C);
  IIC = trace(C * C);
  IIStarC = 0.5 * (IC^2 - IIC);
  IS = trace(S);
  IIS = trace(S' * S);
  IIStarS = 0.5 * (IS^2 - IIS);
 
  % here's symmetric dirichlet 
  dIIStarC = 2 * IC * F - 2 * F * F' * F;
  P = 2 * F + dIIStarC / J^2 - (2 / J^3) * IIStarC * djdf;

  % here we add the ARAP part
  P = P - 2 * ((1 + IS / J) * R - (IIStarS / J^2) * djdf - (1 / J) * F);
  
  P = P/2;
end
