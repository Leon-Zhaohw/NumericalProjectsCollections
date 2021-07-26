%AMG_USERSET_FD sets us up for User Specified Example using FInite
%Difference Discretization

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

amg_globals;

filename = input('\nPlease enter the file name containing your FD matrix.\n(Press Return for the Default):', 's');
if isempty(filename)
    filename = 'anisox_81.txt';
end

fin = fopen (filename,'rt');
dim= fscanf(fin,'%4d',2);
FINEPOINTS = dim(1);
A(1).matrix = fscanf(fin,'%g',[dim(1),dim(2)]); 
fclose(fin);

X_Guess = zeros( dim(2), 1 );

RHS = ones( dim(2), 1 );