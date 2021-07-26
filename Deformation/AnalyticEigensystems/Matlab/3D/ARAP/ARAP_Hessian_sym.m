function [H] = ARAP_Hessian_sym(F)
  % get the twist modes
  T0 = sym([0 -1 0; 1 0 0; 0 0 0]);
  T0 = (sym(1) / sqrt(sym(2))) * T0;
  
  T1 = sym([0 0 0; 0 0 1; 0 -1 0]);
  T1 = (sym(1) / sqrt(sym(2))) * T1;

  T2 = sym([0 0 1; 0 0 0; -1 0 0]);
  T2 = (sym(1) / sqrt(sym(2))) * T2;

  % get the flattened versions
  t0 = vec(T0);
  t1 = vec(T1);
  t2 = vec(T2);

  % get the singular values in an order that is consistent
  % with the numbering in the paper
  s0 = F(1,1);
  s1 = F(2,2);
  s2 = F(3,3);
  
  H = 2 * sym(eye(9,9));
  H = H - (4 / (s0 + s1)) * (t0 * t0');
  H = H - (4 / (s1 + s2)) * (t1 * t1');
  H = H - (4 / (s0 + s2)) * (t2 * t2');
end
