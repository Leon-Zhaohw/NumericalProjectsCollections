function writeOBJ(filename, mesh)
% writes an OBJ file with vertex/face information
%
% Usage:
%   writeOBJ(filename, mesh)

%
% See corresponding paper: "Mixed finite elements for variational surface
% modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
% 2010
%
% (C) 2007 Denis Kovacs, NYU
%-------------------------------------------------------------------------

disp(['writing: ',filename]);
f = fopen( filename, 'w' );

for k=1:size(mesh.V,1)
    fprintf( f, 'v %f %f %f\n', mesh.V(k,1), mesh.V(k,2), mesh.V(k,3) );
end

if isfield(mesh, 'UV')
    for k=1:size(mesh.UV,1)
        fprintf( f, 'vt %f %f\n', mesh.UV(k,1), mesh.UV(k,2) );
    end
end

if isfield(mesh, 'N')
    for k=1:size(mesh.V,1)
        fprintf( f, 'vn %f %f %f\n', mesh.N(k,1), mesh.N(k,2), mesh.N(k,3) );
    end
end
    
    
if ( (~isfield(mesh,'N')) && (~isfield(mesh,'UV')) )
    for k=1:size(mesh.F,1)
        fprintf( f, 'f %d %d %d\n', ...
            mesh.F(k,1), mesh.F(k,2), mesh.F(k,3));
    end
elseif ( (~isfield(mesh,'N')) && (isfield(mesh,'UV')) )
    for k=1:size(mesh.F,1)
        fprintf( f, 'f %d/%d %d/%d %d/%d\n', ...
            mesh.F(k,1), mesh.F(k,1), ...
            mesh.F(k,2), mesh.F(k,2), ...
            mesh.F(k,3), mesh.F(k,3) );
    end
elseif ( (isfield(mesh,'N')) && (~isfield(mesh,'UV')) )
    for k=1:size(mesh.F,1)
        fprintf( f, 'f %d//%d %d//%d %d//%d\n', ...
            mesh.F(k,1), mesh.F(k,1), ...
            mesh.F(k,2), mesh.F(k,2), ...
            mesh.F(k,3), mesh.F(k,3) );
    end
elseif ( (isfield(mesh,'N')) && (isfield(mesh,'UV')) )
    for k=1:size(mesh.F,1)
        fprintf( f, 'f %d/%d/%d %d/%d/%d %d/%d/%d\n', ...
            mesh.F(k,1), mesh.F(k,1), mesh.F(k,1),...
            mesh.F(k,2), mesh.F(k,2), mesh.F(k,2),...
            mesh.F(k,3), mesh.F(k,3), mesh.F(k,3) );
    end
end


fclose(f);
