function [spec_rad] = spectral_radius(level)

% Divya Bansal
% Department of Computer Science
% University of Kentucky

eigenvalues = eig(A(level).matrix);
i =1;
len = size(eigenvalues);
while (i<= len)
    if (eigenvalues(i)<0)
        eigenvalues(i)= eigenvalues(i) * -1;
    end
end

spec_rad = max(eigenvalues);