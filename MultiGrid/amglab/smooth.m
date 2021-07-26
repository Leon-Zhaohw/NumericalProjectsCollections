function [X_apx, iterations] = smooth(level, b, u_in)

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

amg_globals;

if DEBUG == 1
    disp('Starting smooth')
end

if SMOOTHER == JACOBI
    [X_apx, iterations] = jacobi_smoother(level, b, u_in);
elseif SMOOTHER == SOR
    [X_apx, iterations] = SOR_smoother(level, b, u_in);
elseif SMOOTHER == RICHARDSON
    [X_apx, iterations] = richardson_smoother(level, b, u_in);
else
    X_apx = u_in;
    iterations = 0;
end
    

if DEBUG == 1
    disp('Exit smooth')
end