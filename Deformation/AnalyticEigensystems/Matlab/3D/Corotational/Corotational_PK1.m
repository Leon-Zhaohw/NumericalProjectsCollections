function [P] = Corotational_PK1(F)
  [U Sigma V] = svd_rv(F);
  S = V * Sigma * V';
  R = U * V';
  IS = trace(S);
  
  P = 2 * (F - R) + (IS - 3) * R;
end