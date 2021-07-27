fid = fopen('##DIMSFILENAME##');
rows = fread(fid, 1, 'int')
cols = fread(fid, 1, 'int')

fclose(fid);

fid = fopen('##INPUTFILENAME##');
A = fread(fid, [cols rows], 'double');
fclose(fid);

size(A)

tic;
fprintf('Starting SVD\n');
[U, S, V] = svd(A, 'econ');
seconds = toc
principal = (diag(S) .* diag(S));

fid = fopen('##UFILENAME##', 'w');
fwrite(fid, cols, 'int');
fwrite(fid, rows, 'int');
U = U';
fwrite(fid, U, 'double');
fclose(fid);

fid = fopen('##SFILENAME##', 'w');
totalS = size(principal);
fwrite(fid, totalS(1), 'int');
fwrite(fid, principal, 'double');
fclose(fid);

exit;
