% mat_split.m
% Splits a m x n matrix into a n x n x k tensor. Each n x n x 1 tensor
% corresponds to one n x n block of the input matrix. m must be divisble by
% k.
%
% Inputs:
%   A: The matrix to split. [m x n matrix]
%   k: The number of splits to perform. [Integer, m/k must be an integer]
%
% Outputs:
%   split: An n x n x k tensor containing the split matrix.
%
% Author: Aidan Gerkis
% Date: 27/03/2024

function split = mat_split(A, k)
    if mod(length(A(:, 1)), k) ~= 0
        error("Unexpected dimension. The number of rows in matrix A must be divisible by k.");
    end

    n = length(A(:, 1))/k;
    
    split = A(1:n, :); % First split

    for i=1:(k-1) % Get k splices, concatenate to previous slices
        split = cat(3, split, A((i*n + 1):(i + 1)*n, :));
    end
end