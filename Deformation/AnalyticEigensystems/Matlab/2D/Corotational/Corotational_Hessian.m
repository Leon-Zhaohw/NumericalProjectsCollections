function [H] = Corotational_Hessian(F)
  [U Sigma V] = svd_rv(F);

  % get the twist modes
  T = [0 -1; 1 0];
  T = (1 / sqrt(2)) * U * T * V';
  
  % get the flattened versions
  t = vec(T);
  
  % get the rotation mode
  R = (1 / sqrt(2)) * U * V';
  r = vec(R);

  % get the singular values in an order that is consistent
  % with the numbering in the paper
  sigma0 = Sigma(1,1);
  sigma1 = Sigma(2,2);
  
  S = V * Sigma * V';
  IS = trace(S);
  kTerm = (IS - 2);
  
  H = 2 * eye(4,4);
  H = H + 2 * ((kTerm - 2) / (sigma0 + sigma1)) * (t * t');
  H = H + 2 * (r * r');
end
