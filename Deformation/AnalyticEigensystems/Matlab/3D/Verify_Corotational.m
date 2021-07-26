isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave)
    pkg load symbolic;
end
addpath('./util');
addpath('./Corotational');
materialName = sprintf('3D Corotational');

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

syms sigma0 sigma1 sigma2 real;
syms I1 I2 I3 real;
d = sym(3);

Fsym = diag([sigma0 sigma1 sigma2]);
H = Corotational_Hessian_sym(Fsym);

% build the eigenvector gallery
eigenvector_gallery_sym;

% we can characterize the directions orthogonal to rotation
% as a 'pinch' and 'bulge', but these two directions can
% actually be chosen arbitrarily.
p = sym([0 0 0 0 1 0 0 0 -1])';
p = p * sym(1) / sqrt(sym(2));
b = sym([1 0 0 0 -1/sym(2) 0 0 0 -1/sym(2)])';
b = b * sym(1) / sqrt(sym(3) / sym(2));

% verify the eigenpairs
IS = trace(Fsym);
lambdas = zeros(9,1);
lambda(1) = 2 + 2 * ((IS - 3) - 2) / (sigma0 + sigma1);
lambda(2) = 2 + 2 * ((IS - 3) - 2) / (sigma1 + sigma2);
lambda(3) = 2 + 2 * ((IS - 3) - 2) / (sigma0 + sigma2);
lambda(4) = 2 + 3;
lambda(5) = 2;
lambda(6) = 2;
lambda(7) = 2;
lambda(8) = 2;
lambda(9) = 2;
Q = [t0 t1 t2 r l0 l1 l2 p b];

for i = 1:length(lambdas)
    Verify_Eigenpair(H, Q(:,i), lambda(i));
end
fprintf('\n\n');
