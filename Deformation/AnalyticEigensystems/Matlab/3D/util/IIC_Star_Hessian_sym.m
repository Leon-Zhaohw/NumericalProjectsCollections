function [H] = IIC_Star_Hessian_sym(F)
  IC = trace(F' * F);
  H = 2 * IC * sym(eye(9,9));
  f = vec(F);
  H = H + 4 * f * f';

  IIC_H = IIC_Hessian_sym(F);
  H = H - 2 * (IIC_H / 4);
end