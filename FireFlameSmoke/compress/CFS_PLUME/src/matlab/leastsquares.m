fid = fopen('LS.matrix');
rows = fread(fid, 1, 'int');
cols = fread(fid, 1, 'int');
A = fread(fid, [cols rows], 'double');
A = A';
fclose(fid);

fid = fopen('LS.vector');
totalSize = fread(fid, 1, 'int');
b = fread(fid, [totalSize 1], 'double');
fclose(fid);

% do the least squares solve
x = A\b;

fid = fopen('LS_result.vector', 'w');
totalSize = size(x);
fwrite(fid, totalSize(1), 'int');
fwrite(fid, x, 'double');
fclose(fid);

exit;
