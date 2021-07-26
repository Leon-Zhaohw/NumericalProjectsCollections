isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave)
    pkg load symbolic;
end
addpath('./util');
addpath('./Corotational');
materialName = sprintf('2D Corotational');

fprintf('============================================\n');
fprintf('Running numerical tests for %s\n', materialName);
fprintf('============================================\n');
fprintf('Numerically verifying expression for PK1:\n');
Verify_PK1(@Corotational_Psi, @Corotational_PK1);
fprintf('Numerically verifying expression for Hessian:\n');
Verify_Hessian(@Corotational_PK1, @Corotational_Hessian);

fprintf('============================================\n');
fprintf('Running symbolic tests for %s\n', materialName);
fprintf('============================================\n');
syms sigma0 sigma1 real;
Fsym = diag([sigma0 sigma1]);
H = Corotational_Hessian_sym(Fsym);

% build the eigenvector gallery
eigenvector_gallery_sym;

% verify the eigenpairs
IS = trace(Fsym);
lambdas = zeros(4,1);
lambda(1) = 2 + 2 * ((IS - 2) - 2) / (sigma0 + sigma1);
lambda(2) = 2 + 2;
lambda(3) = 2;
lambda(4) = 2;
Q = [t r p l];

for i = 1:length(lambdas)
    Verify_Eigenpair(H, Q(:,i), lambda(i));
end
fprintf('\n\n');
