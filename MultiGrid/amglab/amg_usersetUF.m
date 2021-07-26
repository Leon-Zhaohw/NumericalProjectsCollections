% AMG_USERSET_MAT sets us up for User Specified Matrix

amg_globals;

name = input('\nPlease enter the UF matrix number.\n:');

Problem = UFget(name);

FINEPOINTS = length(Problem(1).A);
A(1).matrix = Problem(1).A;

X_Guess = zeros( FINEPOINTS, 1 );

RHS = ones( FINEPOINTS, 1 );
