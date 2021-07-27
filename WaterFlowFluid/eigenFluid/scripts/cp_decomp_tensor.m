fname = '../Type6Matlab';
in = fopen(fname, 'r');
rows = fread(in, 1, 'int');
cols = fread(in, 1, 'int');
Tensor = zeros(rows,rows,rows);
for i = 1:rows - 1
  Tensor(:,:,i) =  fread(in, [rows, rows], 'double');
  temp1 = fread(in, 1, 'int');
  temp2 = fread(in, 1, 'int');
end
Tensor(:,:, rows - 1) = fread(in, [rows, rows], 'double');
Tensor = tensor(Tensor);
P = parafac_als(Tensor,1000);

