% AMG_GLOBALS contains global variables and flags for the AMG algorithm

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

% Global definitions for yes and no
global YES; YES=1;
global NO; NO = 0;

% Global flags for setup options
global PROB_SRC % What is the source for the problem
    global USER_SPEC_FD; USER_SPEC_FD = 11; % user defined finite difference
    global USER_SPEC_FEM; USER_SPEC_FEM = 12; % user defined finite element
    global EXAMPLE1; EXAMPLE1 = 13; % finite difference poisson
    global USER_SPEC_UF; USER_SPEC_UF = 14; % user defined matrix
    
global PLOT_ACTUAL; % whether or not to calculate and plot actual solution
global SHOW_PROFILE; % whether or not to show the function profiler
    
global PROB_TYPE % How many data points to solve over
    global STIFFNESS; STIFFNESS = 101;
    global ELEMENTS; ELEMENTS = 102;

global SETUP_ALG % Which coarse point selector and intergrid weight generator to use
    global RUGE_STUEBEN; RUGE_STUEBEN = 201;
    global AMGm; AMGm = 202;
    global SMOOTH_AGGREGATE; SMOOTH_AGGREGATE = 203;
    global AMGe; AMGe = 204;
    global BECK; BECK = 205;
    
global CYCLE_TYPE % What "shape" of recursive cycle to use
    global VCYCLE; VCYCLE = 301;
    global WCYCLE; WCYCLE = 302;
    global FCYCLE; FCYCLE = 303;
    
global SETUP_OPT % When to perform the setup
    global AT_ONCE; AT_ONCE = 401;
    global AT_EACH; AT_EACH = 402;

% Global flags for AMG options
global SMOOTHER % which smoother to use
    global JACOBI; JACOBI = 501;
    global SOR; SOR = 502;
    global RICHARDSON; RICHARDSON = 503;

global C_SOLVER % which solver to use on coarsest level
    global THE_SMOOTHER; THE_SMOOTHER = 601;
    global DIRECT_ELIM; DIRECT_ELIM = 602;
    
global POST_SMOOTH % do we post-smooth at each level

global STOP_TYPE % which type of stopping criteria to check for
    global RESID_THRESHOLD; RESID_THRESHOLD = 701;
    global RESID_REDUCE; RESID_REDUCE = 702;
    global MAX_ITRS; MAX_ITRS = 703;
    
global STOP_VALUE; % value for stopping criteria
global ABSOLUTE_MAX_ITRS; % max smoother iterations to execute if stopping criteria not met

% Global variables for use in AMG algorithm
global FINEPOINTS; % number of data points in finest grid (FD or FEM)
global MAXVERTEX; % number of vertices attached to each finite element
global PAIRS; % local node numberings for FEM
global COARSEST; % coarsest level 
global CYCLES; % number of complete recursive cycles to perform

global A; % structure to contain system matrices
% structure elements
% A( x ).matrix - system coefficient matrix at level x (FD)
% A( x ).elements - element matrices for nodes at level x (FEM)
% A( x ).elconn - element connectivity matrix for nodes at level x (FEM)
% A( x ).edges - edge matrices for nodes at level x (FEM)
% A( x ).edconn - edge connectivity matrix for nodes at level x (FEM)

global W; % structure to contain interpolation weights
% structure elements
%  W( x ).Iweight - The interpolation operator between level x+1 (coarse) and x (fine)
%  W( x ).Rweight - The restriction operator between levles x (fine) and x+1 (coarse)
% *NOTE: In most cases, the Iweight and Rweight are transposes of one another
%  W( x ).molecules - The intergrid transfer molecules between levels x and x+1

global SW_BOUND; % boundary between what are considered strong and weak connections
% this factor must be bounded between 0 and 1

global X_Guess; % input initial guess
global RHS; % input right hand side
global SOLN; % computed result

% Global variables to store interesting analysis data
global INIT_RESID; % stores 2-norm of initial residual
global FINAL_RESID; % stores 2-norm of final residual
global RH1; % 2-D stores a history of the 2-norm residual before correction [ RH1(i,j) = resid1 at cycle i level j ]
global RH2; % 2-D stores a history of the 2-norm residual after correction [ RH2(i,j) = resid2 at cycle i level j ]
global IH1; % 2-D stores a history of the number of pre-smoothing iterations [ IH1(i,j) = #itrs at cycle i level j ]
global IH2; % 2-D stores a history of the number of post-smoothing iterations [ IH2(i,j) = #itrs at cycle i level j ]
% Added by Derrick Cerwinsky
global IH;  % 1-D stores a history of the number of smoothing iterations [ IH(j) = #itrs at level j ]
global SHOW_IH;  % This is a flag to show number of smoothing iterations at each level.


% Global variable for storing the SOR relaxation scalar.  Added by Derrick
% Cerwinsky
global OMEGA;

% Global variable for storing the Ruge_Stueben Strong Connection constant.  
% This should range from 0 to 1.   Added by Derrick Cerwinsky
global THETA;

% Global variable for debugging.  Set to 1 for debugging.
global DEBUG;