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

% do the non-negative least squares solve
%x = lsqnonneg(A,b);
[rows cols] = size(A);
x0 = ones(cols,1);
out = bbnnls(A,b,x0);

fid = fopen('LS_result.vector', 'w');
totalSize = size(out.x);
fwrite(fid, totalSize(1), 'int');
fwrite(fid, out.x, 'double');
fclose(fid);

exit;
