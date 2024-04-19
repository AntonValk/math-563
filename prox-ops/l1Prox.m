function P = l1Prox(x, l)
    % l1Prox.m
    %
    % Computes the proximal operator of the l1-norm, as calculated in [1].
    % Accepts as a parameter lambda, the scaling factor of the l2-norm in the
    % Moreau Envelope.
    %
    % Inputs:
    %   x: A matrix in R^(n x n) with which to compute the prox operator. 
    %   l: The scaling factor of the moreau envelope. [OPTIONAL]
    %
    % Outputs:
    %   P: The proximal operator of the l1-norm evaluated at x. [Matrix]
    %
    % Usage:
    %   y = l1Prox(x) outputs the proximal operator of the l1-norm with lambda = 1
    %   y = l1Prox(x, c) outputs the proximal operator of the l1-norm with lambda = c
    %
    % Author: Aidan Gerkis
    % Date: 14-03-2024
    %
    % References:
    %   [1]: C. Paquette, MATH 463/563 Assignment 4 (2024).
    
    % Process inputs
    [m, n] = size(x);

    if nargin == 1 % Process input, determine if a scaling factor is input
        l = 1; % If no scaling factor is input then set l to 1
    end

    % Validate input
    if l < 0 % Check sign of scaling factor
        error("Scaling factor must be positive and non-zero.");
    end

    % Compute proximal operator
    P = sign(x).*max(cat(3, zeros(m, n), abs(x) - l), [], 3);
end