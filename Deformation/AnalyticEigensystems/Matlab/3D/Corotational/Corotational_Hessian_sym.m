function [H] = Corotational_Hessian_sym(F)
  % get the twist modes
  T0 = sym([0 -1 0; 1 0 0; 0 0 0]);
  T0 = (sym(1) / sqrt(sym(2))) * T0;
  
  T1 = sym([0 0 0; 0 0 1; 0 -1 0]);
  T1 = (sym(1) / sqrt(sym(2))) * T1;

  T2 = sym([0 0 1; 0 0 0; -1 0 0]);
  T2 = (sym(1) / sqrt(sym(2))) * T2;

  R = sym([1 0 0; 0 1 0; 0 0 1]);
  R = (sym(1) / sqrt(sym(3))) * R;

  % get the flattened versions
  t0 = vec(T0);
  t1 = vec(T1);
  t2 = vec(T2);
  r = vec(R);

  % get the singular values in an order that is consistent
  % with the numbering in the paper
  sigma0 = F(1,1);
  sigma1 = F(2,2);
  sigma2 = F(3,3);
  
  IS = trace(F);
  kTerm = (IS - 3);
  
  H = 2 * sym(eye(9,9));
  H = H + 2 * ((kTerm - 2) / (sigma0 + sigma1)) * (t0 * t0');
  H = H + 2 * ((kTerm - 2) / (sigma1 + sigma2)) * (t1 * t1');
  H = H + 2 * ((kTerm - 2) / (sigma0 + sigma2)) * (t2 * t2');
  H = H + 3 * (r * r');
end
