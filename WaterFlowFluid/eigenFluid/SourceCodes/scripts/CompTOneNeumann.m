% Numerically compute the tensor for one Neumann wall case.

num_basis_root = 4;
num_basis = num_basis_root * num_basis_root;

C = zeros(16,16,16);

total = 0;

for k = 0:(num_basis - 1)
    for j = 0:(num_basis - 1)
        for i = 0:(num_basis - 1)
            [k1, k2] = Lookup(k, num_basis_root);
            [j1, j2] = Lookup(j, num_basis_root);
            [i1, i2] = Lookup(i, num_basis_root);
            entry = ComputeTensorEntry(i1 - 0.5,i2,j1 - 0.5,j2,k1 - 0.5,k2);
            % Matlab begin with 1.
            C(i+1,j+1,k+1) = entry;
            if (abs(entry) > 1e-10)
                total = total + 1;
            end
        end
    end
end

disp(sprintf('Total number of non-zero entries: %d ', total));

fd = fopen('tensor.dat','w');
fwrite(fd,C,'double');
fclose(fd);