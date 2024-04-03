% imopt_l2.m
% 
% Explicitly computes the L2-norm at a point x. Treats x as a vector.
%
% Inputs:
%   x: The point at which to compute the L2 norm. [m x n Matrix]
%
% Outputs:
%   f: ||x||_2. [Double]
%
% Author: Aidan Gerkis
% Date: 03-04-2024

function f = imopt_l2(x)
    f = sqrt(sum(x.^2, 'all'));
end