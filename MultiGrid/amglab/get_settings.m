function outstring = get_settings()
% Get Settings analyzes the values of the global AMG settings variables,
% writes appropriate descriptions of their state to a string and returns
% that string as "outstring"

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

amg_globals;

% == Problem Information ============================================

p  =[sprintf('\nProblem                                  ')];
if (PROB_SRC==USER_SPEC_FD)
p=[p,sprintf('\n    User Specified FD                    ')];
elseif (PROB_SRC==USER_SPEC_FEM)
p=[p,sprintf('\n    User Specified FEM                   ')];
elseif (PROB_SRC==USER_SPEC_UF)
p=[p,sprintf('\n    User Specified UF Matrix             ')];
elseif (PROB_SRC==EXAMPLE1)
p=[p,sprintf('\n    Example 1 (Poisson)                  ')];
if(PROB_TYPE==STIFFNESS)
p=[p,sprintf('\n    Problem size = %3d x %3d Stiffness   ',FINEPOINTS,FINEPOINTS)];
elseif(PROB_TYPE==ELEMENTS)
p=[p,sprintf('\n    %3d Finite Elements                  ',FINEPOINTS)];
end
else
p=[p,sprintf('\n    Invalid Problem Spec                 ')];
end

% == AMG Setup Algorithm ==============================================

su =  [sprintf('\nAMG Setup                                   ')];
if (SETUP_ALG==RUGE_STUEBEN)
su=[su,sprintf('\n    Ruge Stueben Algorithm                  ')];
elseif (SETUP_ALG==AMGm)
su=[su,sprintf('\n    AMGm Algorithm                          ')];
elseif (SETUP_ALG==SMOOTH_AGGREGATE)
su=[su,sprintf('\n    Smoothed Aggregation Algorithm          ')];
elseif (SETUP_ALG==AMGe)
su=[su,sprintf('\n    AMGe Algorithm                          ')];
elseif (SETUP_ALG==BECK)
su=[su,sprintf('\n    Beck Algorithm                          ')];
else
su=[su,sprintf('\n    Invalid Setup Algorithm                 ')];
end
if (SETUP_OPT==AT_ONCE)
su=[su,sprintf('\n    All Coarse Grids Calculated First       ')];
elseif (SETUP_OPT==AT_EACH)
su=[su,sprintf('\n    Coarse Grids Calculated at Each Level   ')];
else
su=[su,sprintf('\n    Invalid Setup Option                    ')];
end

% == Solving Information ==============================================

s  =[sprintf('\nSolving Information                         ')];
if (SMOOTHER==JACOBI)
s=[s,sprintf('\n    Smoother: Jacobi                        ')];
elseif (SMOOTHER==SOR)
s=[s,sprintf('\n    Smoother: SOR                           ')];
elseif (SMOOTHER==RICHARDSON)
s=[s,sprintf('\n    Smoother: Richardson                    ')];
else
s=[s,sprintf('\n    Invalid Smoother                        ')];
end
if (STOP_TYPE==RESID_THRESHOLD)
s=[s,sprintf('\n    Stop Criteria: Residual Threshold %e    ',STOP_VALUE)];
elseif (STOP_TYPE==RESID_REDUCE)
s=[s,sprintf('\n    Stop Criteria: %e times Initial Residual',STOP_VALUE)];
elseif (STOP_TYPE==MAX_ITRS)
s=[s,sprintf('\n    Stop Criteria: %3d Iterations           ',STOP_VALUE)];
else
s=[s,sprintf('\n    Invalid Stopping Criteria               ')];
end
s=[s,sprintf('\n    Absolute Iteration Limit: %3d           ',ABSOLUTE_MAX_ITRS)];
if (C_SOLVER==THE_SMOOTHER)
s=[s,sprintf('\n    Coarse Level Solver: Current Smoother   ')];
elseif (C_SOLVER==DIRECT_ELIM)
s=[s,sprintf('\n    Coarse Level Solver: Direct Elimination ')];
%elseif (C_SOLVER==LU_DECOMP)
%s=[s,sprintf('\n    Coarse Level Solver: LU Decomposition   ')];
else
s=[s,sprintf('\n    Invalid Coarse Level Solver             ')];    
end
if (POST_SMOOTH==YES)
s=[s,sprintf('\n    Post-Smoothing Activated                ')];
end

% == Multigrid Settings ===========================================

m  =[sprintf('\nMultigrid Settings                          ')];

if (CYCLE_TYPE==VCYCLE)
m=[m,sprintf('\n    Cycle Type: V-Cycle                     ')];
elseif (CYCLE_TYPE==VCYCLE)
m=[m,sprintf('\n    Cycle Type: W-Cycle                     ')];
elseif (CYCLE_TYPE==VCYCLE)
m=[m,sprintf('\n    Cycle Type: F-Cycle                     ')];
else
m=[m,sprintf('\n    Invalid Cycle Type                      ')];
end
m=[m,sprintf('\n    Number of Levels = %1d                  ',COARSEST)];
m=[m,sprintf('\n    Number of Cycles = %1d                  ',CYCLES)];

outstring = sprintf('%s\n%s\n%s\n%s\n',p,su,s,m);
