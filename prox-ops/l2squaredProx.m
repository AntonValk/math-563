% l2squaredProx.m
%
% Computes the proximal operator of the l2_squared-norm, as calculated in [1].
% Accepts as a parameter lambda, the scaling factor of the l2-norm in the
% Moreau Envelope.
%
% Inputs:
%   x: A matrix in R^(n x n) with which to compute the prox operator. 
%   l: The scaling factor of the moreau envelope. [OPTIONAL]
%
% Outputs:
%   P: The proximal operator of the l2_squared-norm evaluated at x. [Matrix]
%
% Usage:
%   y = l1Prox(x) outputs the proximal operator of the l1-norm with lambda = 1
%   y = l1Prox(x, c) outputs the proximal operator of the l1-norm with lambda = c
%
% Author: Linda Hu
% Date: 21-03-2024
%
% References:
%   [1]: C. Paquette, MATH 463/563 Assignment 4 (2024).
%   [2]: Aiden Gerkis' code for l1prox.m

function P = l2squaredProx(x, l)
    switch nargin % Process input, determine if a scaling factor is input
        case 1 % If no scaling factor is input then set l to 1
            l = 1;
        otherwise
            l = l;
    end

    % Validate input
    if l < 0 % Check sign of scaling factor
        error("Scaling factor must be positive.");
    end

    % Compute proximal operator
    P = x/(2*l+1);
end