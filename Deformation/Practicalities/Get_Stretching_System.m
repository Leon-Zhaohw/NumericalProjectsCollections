function [A] = Get_Stretching_System(Psi)
  syms I1 I2 I3 sigma0 sigma1 sigma2 sigmai sigmak real;
  A = sym(zeros(3,3));

  % first derivatives
  firstI1 = diff(Psi, I1);
  firstI2 = diff(Psi, I2);
  firstI3 = diff(Psi, I3);

  % second derivatives
  secondI1 = diff(firstI1, I1);
  secondI2 = diff(firstI2, I2);
  secondI3 = diff(firstI3, I3);

  % mixed derivatives
  secondI1I2 = diff(firstI1, I2);
  secondI1I3 = diff(firstI1, I3);
  secondI2I3 = diff(firstI2, I3);

  % get the diagonal entry
  aii = 2 * firstI2 + secondI1 + 4 * sigmai^2 * secondI2 + ...
        (I3 / sigmai)^2 * secondI3 + 4 * sigmai * secondI1I2 + ...
        4 * I3 * secondI2I3 + 2 * (I3 / sigmai) * secondI1I3;
  A(1,1) = subs(aii, sigmai, sigma0);
  A(2,2) = subs(aii, sigmai, sigma1);
  A(3,3) = subs(aii, sigmai, sigma2);

  % get the off-diagonal entry
  aij = sigmak * firstI3 + secondI1 + 4 * (I3 / sigmak) * secondI2 + ...
        sigmak * I3 * secondI3 + ...
        2 * sigmak * (I2 - sigmak^2) * secondI2I3 + ...
        (I1 - sigmak) * (sigmak * secondI1I3 + 2 * secondI1I2);
  A(1,2) = subs(aij, sigmak, sigma2);
  A(1,3) = subs(aij, sigmak, sigma1);
  A(2,3) = subs(aij, sigmak, sigma0);

  % symmetrize and simplify
  A(2,1) = A(1,2);
  A(3,1) = A(1,3);
  A(3,2) = A(2,3);
  A = Simplify_Invariants(A);
end
