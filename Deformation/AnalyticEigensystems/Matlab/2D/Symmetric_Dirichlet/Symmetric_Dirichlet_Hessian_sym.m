function [H] = Symmetric_Dirichlet_Hessian_sym(F)
  djdf = DJDF(F);
  g = vec(djdf);
  f = vec(F);

  J = det(F);
  J2 = J * J;

  anti = sym(fliplr(diag([1 -1 -1 1])));
  I2 = trace(F' * F);

  H = (1 + 1 / J2) * eye(4,4);
  H = H - (I2 / J^3) * anti;
  H = H - (2 / J^3) * (g * f' + f * g');
  H = H + (3 * I2 / J^4) * (g * g');
end
