function [H] = IIC_Hessian_sym(F)
  H = sym(zeros(9,9));
  for i = 1:9
    [DF] = DFDF(i);
    A = 4 * (DF * F' * F + F * F' * DF + F * DF' * F);
    column = reshape(A, 9,1);
    H(:,i) = column;
  end
end
