function mesh = readOFF( filename )
% reads an OFF file with vertex/face information
%
% Usage:
%   mesh = readOFF(filename)

%
% See corresponding paper: "Mixed finite elements for variational surface
% modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
% 2010
%
% (C) 2007 Denis Kovacs, NYU
%-------------------------------------------------------------------------

mesh = [];
fp = fopen( filename, 'r' );
OFFheader = upper(fscanf( fp, '%s\n', 1 ));
if (OFFheader(end-2:end) ~= 'OFF') warning('no OFF file!'); return; end
OFFdim = 3;
OFF_N = 0; OFF_C=0; OFF_ST=0;

if find(OFFheader=='N') OFFdim = OFFdim+3; OFF_N=1; end
if find(OFFheader=='C') OFFdim = OFFdim+3; OFF_C=1; end
if find(OFFheader=='S') OFFdim = OFFdim+2; OFF_ST=1; end

d = fscanf( fp, '%d', 3);
nV = d(1); nF = d(2); nE = d(3);

disp(sprintf('  - Reading %d vertices', nV));

switch OFFdim
    case  3; OFFV = textscan( fp, '%f %f %f', nV);
    case  5; OFFV = textscan( fp, '%f %f %f %f %f', nV);
    case  6; OFFV = textscan( fp, '%f %f %f %f %f %f', nV);
    case  7; OFFV = textscan( fp, '%f %f %f %f %f %f %f', nV);
    case  8; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f', nV);
    case  9; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f %f', nV);
    case 10; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f %f %f', nV);
    case 11; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f %f %f %f', nV);
    otherwise; error('Unsupported number of vertex entries');
end

try
   OFFV = cell2mat(OFFV); 
end

OFFdim = 1;
mesh.V = OFFV(:,OFFdim:(OFFdim+2)); OFFdim = OFFdim + 3;
if (OFF_N) mesh.N = OFFV(:,OFFdim:(OFFdim+2)); OFFdim = OFFdim + 3; end
if (OFF_C) mesh.C = OFFV(:,OFFdim:(OFFdim+2)); OFFdim = OFFdim + 3; end
if (OFF_ST) mesh.UV = OFFV(:,OFFdim:(OFFdim+1)); OFFdim = OFFdim + 2; end


disp(sprintf('  - Reading %d faces', nF));
F = cell2mat( textscan( fp, '%d %d %d %d %d %d', nF ) );
mesh.F = double( F(:, 2:(F(1,1)+1) ) + 1 );

disp('  - done.');
