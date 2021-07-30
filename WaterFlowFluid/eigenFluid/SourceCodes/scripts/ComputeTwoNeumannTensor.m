num_basis_root = 5;
num_basis = num_basis_root * num_basis_root;

C = zeros(25,25,25);

total = 0;

for k = 0:(num_basis - 1)
    for j = 0:(num_basis - 1)
        for i = 0:(num_basis - 1)
            [k1, k2] = Lookup(k, num_basis_root);
            [j1, j2] = Lookup(j, num_basis_root);
            [i1, i2] = Lookup(i, num_basis_root);
            entry = ComputeTwoNeumannTensorEntry(i1 - 1,i2 - 1,j1 - 1,j2 - 1,k1 - 1,k2 - 1);
            % Matlab begin with 1.
            C(i+1,j+1,k+1) = entry;
            if (abs(entry) > 1e-10)
                total = total + 1;
            end
        end
    end
end

disp(sprintf('Total number of non-zero entries: %d ', total));

fd = fopen('tensor_2N.dat','w');
fwrite(fd,C,'double');
fclose(fd);