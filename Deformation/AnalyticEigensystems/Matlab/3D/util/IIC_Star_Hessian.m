function [H] = IIC_Star_Hessian(F)
  IC = trace(F' * F);
  H = 2 * IC * eye(9,9);
  f = vec(F);
  H = H + 4 * f * f';

  IIC_H = IIC_Hessian(F);
  H = H - 2 * (IIC_H / 4);
end
