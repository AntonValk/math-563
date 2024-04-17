% imopt_psnr.m
%
% Computes the Peak Signal to Noise Ratio (PSNR)
%
% Inputs:
%   x: A matrix in R^(n x n) that represents the ground truth image.
%   x_hat: A matrix in R^(n x n) that represent the reconstructed image.
%
% Outputs:
%   e: Peak Signal to Noise Ratio (PSNR) (scalar)
%
% Usage:
%   e = imopt_psnr(x, x_hat) outputs the PSNR of the reconstructed image with
%   respect to ground truth
%
% Author: Antonios Valkanas
% Date: 02-04-2024
%
% References:
%   [1]: C. Paquette, MATH 463/563 Assignment 4 (2024).


function e = imopt_psnr(x, x_hat)
    [N, M] = size(x);
    mse = 1/(N*M)*sum((x - x_hat).^2, 'all');
    e = -1 * log(mse);
end
