function mesh = readOBJ( filename )
% reads an OBJ file with vertex/face information
%
% Usage:
%   mesh = readOBJ(filename)

%
% See corresponding paper: "Mixed finite elements for variational surface
% modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
% 2010
%
% (C) 2007 Denis Kovacs, NYU
%-------------------------------------------------------------------------

mesh = [];
V = [];
UV = [];
F = {};
fp = fopen( filename, 'r' );
type = fscanf( fp, '%s', 1 );
while strcmp( type, '' ) == 0
    if strcmp( type, 'v' ) == 1
        v = fscanf( fp, '%g %g %g\n' );
        V = [V; v'];
    elseif strcmp( type, 'vt')
        v = fscanf( fp, '%g %g %g\n' );
        UV = [UV; v'];
    elseif strcmp( type, 'f' ) == 1
        line = fgets(fp);
        [t, count] = sscanf(line, '%d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d');

        if (count>2)
            t = t(1:3:end);
        else
            [t, count] = sscanf(line, '%d/%d %d/%d %d/%d %d/%d %d/%d');
            if (count>2)
                t = t(1:2:end);
            else
                [t, count] = sscanf( line, '%d %d %d %d %d %d %d %d %d %d %d\n' );
            end
        end
        F = [F; {t'}];
    elseif strcmp( type, '#' ) == 1
        fscanf( fp, '%s\n', 1 );
    end
    type = fscanf( fp, '%s', 1 );
end
fclose( fp );

try
    F = cell2mat(F);
end

%% transform into array if all faces have the same number of vertices

mesh.V = V;
if (size(UV,1)>0) mesh.UV = UV; end
mesh.F = F;

end
