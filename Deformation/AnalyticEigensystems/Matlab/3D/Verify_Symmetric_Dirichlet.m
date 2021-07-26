isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave)
    pkg load symbolic;
end
addpath('./util');
addpath('./Symmetric_Dirichlet');
materialName = sprintf('3D Symmetric_Dirichlet');

fprintf('============================================\n');
fprintf('Running numerical tests for %s\n', materialName);
fprintf('============================================\n');
fprintf('Numerically verifying expression for PK1:\n');
Verify_PK1(@Symmetric_Dirichlet_Psi, @Symmetric_Dirichlet_PK1);
fprintf('Numerically verifying expression for Hessian:\n');
Verify_Hessian(@Symmetric_Dirichlet_PK1, @Symmetric_Dirichlet_Hessian);

fprintf('============================================\n');
fprintf('Running symbolic tests for %s\n', materialName);
fprintf('============================================\n');

syms sigma0 sigma1 sigma2 real;
syms I1 I2 I3 real;
d = sym(3);

Fsym = diag([sigma0 sigma1 sigma2]);
H = Symmetric_Dirichlet_Hessian_sym(Fsym);

Psi = I2 / 2 + (1 / sym(8)) * ((I1^2 - I2) / I3)^2 - I1 / I3;

[Q lambdas] = Get_Analytic_Eigensystem(Psi);

for i = 1:length(lambdas)
    Verify_Eigenpair(H, Q(:,i), lambdas(i));
end
fprintf('\n\n');
