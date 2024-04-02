% rmse.m
%
% Computes the Root Mean Squared Error (RMSE)
%
% Inputs:
%   x: A matrix in R^(n x n) that represents the ground truth image.
%   x_hat: A matrix in R^(n x n) that represent the reconstructed image.
%
% Outputs:
%   e: Root Mean Squared Error (RMSE) (scalar)
%
% Usage:
%   e = rmse(x, x_hat) outputs the RMSE of the reconstructed image with
%   respect to ground truth
%
% Author: Antonios Valkanas
% Date: 02-04-2024
%
% References:
%   [1]: C. Paquette, MATH 463/563 Assignment 4 (2024).


function e = imopt_rmse(x, x_hat)
%ERRORS Summary of this function goes here
%   Detailed explanation goes here
[N, ~] = size(x);
[~, M] = size(x);
e = 1/(N*M)*sum((x - x_hat).^2, 'all');
end

