function [H] = Corotational_Hessian(F)
  [U Sigma V] = svd_rv(F);

  % get the twist modes
  T0 = [0 -1 0; 1 0 0; 0 0 0];
  T0 = (1 / sqrt(2)) * U * T0 * V';
  
  T1 = [0 0 0; 0 0 1; 0 -1 0];
  T1 = (1 / sqrt(2)) * U * T1 * V';

  T2 = [0 0 1; 0 0 0; -1 0 0];
  T2 = (1 / sqrt(2)) * U * T2 * V';

  % get the flattened versions
  t0 = vec(T0);
  t1 = vec(T1);
  t2 = vec(T2);
  
  % get the rotation mode
  R = (1 / sqrt(3)) * U * V';
  r = vec(R);

  % get the singular values in an order that is consistent
  % with the numbering in the paper
  sigma0 = Sigma(1,1);
  sigma1 = Sigma(2,2);
  sigma2 = Sigma(3,3);
  
  S = V * Sigma * V';
  IS = trace(S);
  kTerm = (IS - 3);
  
  H = 2 * eye(9,9);
  H = H + 2 * ((kTerm - 2) / (sigma0 + sigma1)) * (t0 * t0');
  H = H + 2 * ((kTerm - 2) / (sigma1 + sigma2)) * (t1 * t1');
  H = H + 2 * ((kTerm - 2) / (sigma0 + sigma2)) * (t2 * t2');
  H = H + 3 * (r * r');
end
