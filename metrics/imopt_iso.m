function f = imopt_iso(x)
    % imopt_iso.m
    % 
    % Explicitly computes the iso-norm at a point (x, y):
    %   ||(x, y)||_ISO = sum(sqrt(x_i^2 + y_i^2))
    %
    % Inputs:
    %   x: The point at which to compute the iso norm. [m x n x 2 Tensor]
    %
    % Outputs:
    %   f: ||(x, y)||_ISO. [Double]
    %
    % Author: Aidan Gerkis
    % Date: 03-04-2024
    
    x1 = x(:, :, 1);
    x2 = x(:, :, 2);
    f = sum(sqrt(x1.^2 + x2.^2), 'all');
end