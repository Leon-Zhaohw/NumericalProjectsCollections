function res = residual(level, b, u_apx)

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

amg_globals;

if DEBUG == 1
    disp('Starting residual')
end

% in Au=B
% the residual is B - A*u_apx
% where B is RHS
% A is system matrix
% and u_apx is the approximated solution vector
res = b - (A(level).matrix * u_apx);

if DEBUG == 1
    disp('Exit residual')
end