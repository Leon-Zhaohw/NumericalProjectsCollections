function fine_resid = interpolate(level, u_coarse)

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

amg_globals;

if DEBUG == 1
    disp('Starting interpolate')
end

fine_resid = W(level).Iweight * u_coarse;

if DEBUG == 1
    disp('Exit interpolate')
end