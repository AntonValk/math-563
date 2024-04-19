function f = imopt_inds(x)
    % imopt_inds.m
    % 
    % Explicitly computes the indicator of x on the component-wise set [0,1]
    %   f = 0 if 0 <= x_i <= 1 for all i
    %   f = Inf otherwise
    %
    % Inputs:
    %   x: The point at which to compute the indicator. [m x n Matrix]
    %
    % Outputs:
    %   f: The indicator of s. [Double]
    %
    % Author: Aidan Gerkis
    % Date: 03-04-2024
    
    if max(x, [], 'all') > 1 || min(x, [], 'all') < 0
        f = Inf;
    else % All values are in the desired range
        f = 0;
    end
end