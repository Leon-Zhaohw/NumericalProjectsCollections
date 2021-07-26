function [H] = Symmetric_ARAP_Hessian(F)
  [U Sigma V] = svd_rv(F);

  R = U * V';
  S = V * Sigma * V';

  djdf = DJDF(F);
  g = vec(djdf);
  f = vec(F);
  r = vec(R);

  J = det(F);
  J2 = J * J;
  J3 = J * J * J;

  anti = fliplr(diag([1 -1 -1 1]));
  I1 = trace(S);
  I2 = trace(F' * F);

  invSqrt2 = 1 / sqrt(2);
  T = [0 1; -1 0];
  T = invSqrt2 * T;
  q = vec(U * T * V');

  H = 2 * (1 + 1 / J2) * eye(4,4);
  H = H - (4 / J3) * (g * f' + f * g');
  H = H + (2 / J2) * (g * r' + r * g');
  H = H + (2 / J2) * (I1 - I2 / J) * anti;
  H = H + (2 / J3) * (3 * I2 / J - 2 * I1) * (g * g');
  H = H - (4 / I1) * (1 + 1 / J) * (q * q');
end
