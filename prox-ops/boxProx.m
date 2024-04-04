% boxProx.m
%
% Computes the component-wise projection onto [0, 1], as defined in [1]. The 
% output is a vector y with y_i in [0, 1].
%
% Inputs:
%   x: A matrix in R^(n x n) with which to compute the projection.
%   l: A scaling factor, if used the output is the component-wise
%      projection onto [0, lambda]. [OPTIONAL]
% Outputs:
%   y: The component-wise projection of x onto [0, 1], in R^(n x n). [Matrix]
%      Returns 0 if an error occurs.
%
% Usage:
%   y = boxProx(x) outputs the component-wise projection of x onto [0, 1]
%   y = boxProx(x, c) outputs the component-wise projection of x onto [0, c]
%
% Author: Aidan Gerkis
% Date: 14-03-2024
%
% References:
%   [1]: C. Paquette, MATH 463/563 Final Project Description (2024).

function y = boxProx(x, l)
    % Process inputs
    [m, n] = size(x);

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
    
    % Compute Box Prox

    y = min(cat(3, max(cat(3, zeros(m, n), x), [], 3), l*ones(m, n)), [], 3);
end