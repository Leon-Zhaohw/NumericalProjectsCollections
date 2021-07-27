function [i1, i2] = Lookup( i,num_basis_root )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if (i > num_basis_root * num_basis_root)
    i1 = -1;
    i2 = -1;
else
    i1 = rem(i, num_basis_root) + 1;
    i2 = fix(i/num_basis_root) +1;
end

end

