function X_soln = gauss_elim(level, b)

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

amg_globals;

if DEBUG == 1
    disp('Starting gauss_elim')
end

% Solves the system A(level) x_soln = b
% This command will use Gauss Elimination when appropriate 
% When not appropriate i.e. A(level) is near singular or badly scaled
% this command will use QR Decomposition with Pivoting

X_soln = A(level).matrix \ b;

if DEBUG == 1
    disp('Exit gauss_elim')
end