function [x] = vec(A)
    [rows cols] = size(A);
    x = reshape(A, rows * cols, 1);
end