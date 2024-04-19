function f = imopt_l1(x)
    % imopt_l1.m
    % 
    % Explicitly computes the L1-norm at a point x. Treats x as a vector.
    %
    % Inputs:
    %   x: The point at which to compute the L1 norm. [m x n Matrix]
    %
    % Outputs:
    %   f: ||x||_1. [Double]
    %
    % Author: Aidan Gerkis
    % Date: 03-04-2024
    
    f = sum(abs(x), 'all');
end