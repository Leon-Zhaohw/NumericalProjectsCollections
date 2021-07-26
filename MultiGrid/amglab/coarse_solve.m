function X_soln = coarse_solve(level, b)

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

amg_globals;

if DEBUG == 1
    disp('Starting coarse_solve')
end

initGuess = zeros(size(b));
if C_SOLVER == THE_SMOOTHER
    X_soln = smooth(level, b, initGuess);
elseif C_SOLVER == DIRECT_ELIM
    X_soln = gauss_elim(level, b);
else
    X_soln = initGuess;
end

if DEBUG == 1
    disp('Exit coarse_solve')
end