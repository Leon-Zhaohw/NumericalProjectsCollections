function [H] = Symmetric_ARAP_Hessian_sym(F)
  djdf = DJDF(F);
  g = vec(djdf);
  f = vec(F);
  R = sym(eye(2,2));
  r = vec(R);

  J = det(F);
  J2 = J * J;
  J3 = J * J * J;

  anti = sym(fliplr(diag([1 -1 -1 1])));
  I1 = trace(F);
  I2 = trace(F' * F);
  
  invSqrt2 = sym(1) / sqrt(sym(2));
  T = sym([0 -1; 1 0]);
  T = invSqrt2 * T;
  q = vec(T);

  H = 2 * (1 + 1 / J2) * eye(4,4);
  H = H - (4 / J3) * (g * f' + f * g');
  H = H + (2 / J2) * (g * r' + r * g');
  H = H + (2 / J2) * (I1 - I2 / J) * anti;
  H = H + (2 / J3) * (3 * I2 / J - 2 * I1) * (g * g');
  H = H - (4 / I1) * (1 + 1 / J) * (q * q');
end
