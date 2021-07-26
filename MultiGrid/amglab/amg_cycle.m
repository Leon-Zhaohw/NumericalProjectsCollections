% AMG_CYCLE Algebraic Multigrid Cycle algorithm.
%
%       U_OUT = AMG_CYCLE(CYCLE, LEVEL, B, U_IN) uses the appropriate cycle to recursively 
%       solve the linear system AX=B at the given level. CYCLE gives
%       the number of the current cycle.
%
%   NOTE: only the V-Cycle is implemented here at this time

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

function u_out = amg_vycle(cycle, level, b, u_in)

amg_globals;

if level == COARSEST % if current cycle is at coarsest level
   RH1(cycle,level) = norm(residual(level, b, u_in)); % store the residual before solution in the history
   u_out            = coarse_solve(level, b); % execute coarsest grid solve
   RH2(cycle,level) = norm(residual(level, b, u_out)); % store the final residual in the history
   IH(level) = 0;% IH2(cycle,level) = 0; % there is no smoothing on this level
else % otherwise...
   [u_apx,numitr]   = smooth(level, b, u_in);  %  smooth the problem at the current level
   IH(level) = IH(level) + numitr; % store the number of iterations in the history
   resid            = residual(level, b, u_apx); % calculate the residual before solution
   RH1(cycle,level) = norm(resid); % store this residual in the history
   apx_course       = restrict(level, u_apx); % restrict approximation to next course level
   b_course         = restrict(level, resid); % restrict residual to next course level
   u_course         = amg_cycle(cycle, level+1, b_course, apx_course); % use cycle to solve residual on next course level
   correct          = interpolate(level, u_course); % interpolate result from course cycle to this level
   u_apx            = u_apx + correct; % use interpolated correction factor to correct approximation
   if POST_SMOOTH % if we are post-smoothing at each level
        [u_out,numitr] = smooth(level, b, u_apx); % smooth the problem at current level using corrected approximation
        IH(level) = IH(level) +  numitr; % store the number of iterations in the history
   else u_out = u_apx; % otherwise, copy the corrected approximation into the return variable
   end
   RH2(cycle,level) = norm(residual(level, b, u_out)); % store residual after correction in the history
end
