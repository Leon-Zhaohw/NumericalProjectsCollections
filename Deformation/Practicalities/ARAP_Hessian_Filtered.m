function [H] = ARAP_Hessian_Filtered(F)
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

  % get the singular values in an order that is consistent
  % with the numbering in the paper
  s0 = Sigma(1,1);
  s1 = Sigma(2,2);
  s2 = Sigma(3,3);

  % build the filter the non-trivial eigenvalues 
  lambda0 = 2 / (s0 + s1);
  lambda1 = 2 / (s1 + s2);
  lambda2 = 2 / (s0 + s2);
  if (s0 + s1 < 2)
    lambda0 = 1;
  end
  if (s1 + s2 < 2)
    lambda1 = 1;
  end
  if (s0 + s2 < 2)
    lambda2 = 1;
  end

  H = eye(9,9);
  H = H - lambda0 * (t0 * t0');
  H = H - lambda1 * (t1 * t1');
  H = H - lambda2 * (t2 * t2');

  H = 2 * H;
end
