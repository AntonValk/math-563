% isoProx.m
%
% Computes the proximal operator of the iso-norm, as defined in [1].
% Accepts as a parameter gamma, the scaling factor of the iso-norm.
%
% Inputs:
%   w: Two matrices in R^(n x n) with which to compute the prox operator.
%      w(:,:,1) = x, w(:,:,2) = y.
%   g: The scaling factor of the iso-norm. [OPTIONAL]
%
% Outputs:
%   P: The proximal operator of the iso-norm evaluated at x. [Matrix]
%
% Usage:
%   y = isoProx(x) outputs the proximal operator of the iso-norm with lambda = 1
%   y = isoProx(x, a) outputs the proximal operator of the iso-norm with gamma = a
%
% Author: Aidan Gerkis
% Date: 14-03-2024
%
% References:
%   [1]: C. Paquette, MATH 463/563 Final Project Description (2024).

function P = isoProx(w, g)
    % Process inputs
    n = length(w(:, 1));

    switch nargin % Process input, determine if scaling factors are input
        case 1 % If no scaling factors are input then set g to 1
            g = 1;
        otherwise 
            g = g;
    end

    % Validate input
    if g <= 0 % Check sign of scaling factor
        disp("Error in isoProx: Scaling factor must be positive and non-zero.")
        P = 0;
        return
    end
    
    % Compute prox op
    x = w(:,:,1);
    y = w(:,:,2);

    P = max(cat(3, zeros(n), sign(sqrt(x.^2 + y.^2) - g)), [], 3)*(1 - g./sqrt(x.^2 + y.^2)).*w;
end