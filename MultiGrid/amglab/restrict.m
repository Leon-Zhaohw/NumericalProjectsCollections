function coarse_resid = restrict(level, resid)

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

amg_globals;

if DEBUG == 1
    disp('Starting restrict')
end

if SETUP_OPT==AT_EACH
    amg_setup(level);
end

% generate coarse grid residual
coarse_resid = W(level).Rweight * resid;

if DEBUG == 1
    disp('Exit restrict')
end