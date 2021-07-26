% AMG_DEFAULTS sets up control flags to default values

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

function amg_defaults;

amg_globals;

PROB_SRC = EXAMPLE1;
PLOT_ACTUAL = NO;
SHOW_PROFILE = NO;
PROB_TYPE = STIFFNESS;
SETUP_ALG = RUGE_STUEBEN;
CYCLE_TYPE = VCYCLE;
SW_BOUND = 0.7;
SETUP_OPT = AT_ONCE;
SMOOTHER = SOR;
C_SOLVER = DIRECT_ELIM;
POST_SMOOTH = NO;
STOP_TYPE = RESID_THRESHOLD;
STOP_VALUE = 1.0e-7;
ABSOLUTE_MAX_ITRS = 100;
FINEPOINTS = 9;
COARSEST = 2; 
CYCLES = 1;
OMEGA = 1.0;
SHOW_IH = NO;
THETA = .25;
DEBUG = 0;