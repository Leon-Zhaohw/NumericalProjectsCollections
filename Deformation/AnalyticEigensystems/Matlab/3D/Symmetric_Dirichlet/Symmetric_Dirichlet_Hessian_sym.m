function [H] = Symmetric_Dirichlet_Hessian_sym(F)
  djdf = DJDF_sym(F);
  g = vec(djdf);

  J = det(F);
  C = F' * F;
  IC = trace(C);
  IIC = trace(C * C);

  IIStar = (IC * IC - IIC) / 2;

  dIIStar = 2 * IC * F - 2 * F * F' * F;
  t = vec(dIIStar);
  d2IIStar = IIC_Star_Hessian_sym(F);
  d2J = HessianJ(F);

  H = 2 * eye(9,9);
  H = H - (2 / J^3) * (g * t' + t * g');
  H = H + (6 * IIStar / J^4) * (g * g');
  H = H + (1 / J^2) * d2IIStar;
  H = H - (2 * IIStar / J^3) * d2J;
  
  H = H / 2;
end
