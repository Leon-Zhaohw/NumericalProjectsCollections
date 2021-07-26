function [H] = ARAP_Hessian(F)
  [U Sigma V] = svd_rv(F);

  invSqrt2 = 1 / sqrt(2);
  T = [0 1; -1 0];
  T = invSqrt2 * T;
  t = vec(U * T * V');

  % get the singular values in an order that is consistent
  % with the numbering in the paper
  s0 = Sigma(1,1);
  s1 = Sigma(2,2);
  
  H = 2 * eye(4,4);
  H = H - (4 / (s0 + s1)) * (t * t');
end
