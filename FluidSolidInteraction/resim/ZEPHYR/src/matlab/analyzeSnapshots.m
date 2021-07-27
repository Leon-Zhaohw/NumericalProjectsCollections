fid = fopen('snapshots.matrix');
rows = fread(fid, 1, 'int')
cols = fread(fid, 1, 'int')
A = fread(fid, [cols rows], 'double');
A = A';
fclose(fid);

tic;
[U, S, V] = svd(A, 'econ');
seconds = toc
principal = (diag(S) .* diag(S));

fid = fopen('U.matrix.matlab', 'w');
fwrite(fid, rows, 'int');
fwrite(fid, cols, 'int');
U = U';
fwrite(fid, U, 'double');
fclose(fid);

fid = fopen('S.vector.matlab', 'w');
totalS = size(principal);
fwrite(fid, totalS(1), 'int');
fwrite(fid, principal, 'double');
fclose(fid);

%exit;