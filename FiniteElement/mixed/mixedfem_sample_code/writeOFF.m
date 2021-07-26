function writeOFF(filename, mesh, format)
% writes an OFF file with vertex/face information
%
% Usage:
%   writeOFF(filename, mesh, (format))
%
% Input:
%   format: optional. ['N','C','ST'] depending on which additional vertex attributes should be stored
%

%
% See corresponding paper: "Mixed finite elements for variational surface
% modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
% 2010
%
% (C) 2008 Denis Kovacs / NYU

disp(['writing: ',filename]);

if (nargin<3) 
    format=[];
    if isfield(mesh,'N') format = ['N',format]; end
    if isfield(mesh,'C') format = ['C',format]; end
    if isfield(mesh,'UV') format = ['UV',format]; end
end

OFFheader = 'OFF';

OFFV = mesh.V;

if find(format=='N')
    OFFheader = ['N',OFFheader];
    if ~isfield(mesh,'N') OFFV = [OFFV, zeros(size(mesh.V))]; else OFFV = [OFFV, mesh.N]; end
end

if find(format=='C') 
    OFFheader = ['C',OFFheader];
    if ~isfield(mesh,'C') OFFV = [OFFV, zeros(size(mesh.V))]; else OFFV = [OFFV, mesh.C]; end
end

if find(format=='U') 
    OFFheader = ['ST',OFFheader];
    if ~isfield(mesh,'UV') OFFV = [OFFV, zeros(size(mesh.V,1),2)]; else OFFV = [OFFV, mesh.UV]; end
end 

f = fopen( filename, 'wt' );
fprintf(f, [OFFheader,'\n']);
fprintf(f, '%d %d 0\n', size(mesh.V,1), size(mesh.F, 1));

switch size(OFFV, 2)
    case  3; fprintf(f, '%f %f %f\n', OFFV');
    case  5; fprintf(f, '%f %f %f %f %f\n', OFFV');
    case  6; fprintf(f, '%f %f %f %f %f %f\n', OFFV');
    case  7; fprintf(f, '%f %f %f %f %f %f %f\n', OFFV');
    case  8; fprintf(f, '%f %f %f %f %f %f %f %f\n', OFFV');
    case  9; fprintf(f, '%f %f %f %f %f %f %f %f %f\n', OFFV');
    case 10; fprintf(f, '%f %f %f %f %f %f %f %f %f %f\n', OFFV');
    case 11; fprintf(f, '%f %f %f %f %f %f %f %f %f %f %f\n', OFFV');
    otherwise; error('Unsupported number of vertex entries');
end

if (~isempty(mesh.F)) mesh.F; mesh.F = mesh.F - 1; end

switch size(mesh.F, 2)
    case 0;
    case 1; fprintf( f, '1 %d\n', mesh.F');
    case 2; fprintf( f, '2 %d %d\n', mesh.F');
    case 3; fprintf( f, '3 %d %d %d\n', mesh.F');
    case 4; fprintf( f, '4 %d %d %d %d\n', mesh.F');
    case 5; fprintf( f, '5 %d %d %d %d %d\n', mesh.F');
    otherwise; error('Unsupported number of vertex entries');
end    

fclose(f);
