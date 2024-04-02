% prox_tg.m
%
%
% Inputs:
%   x: The point at which to evaluate the prox. [n x n x 3 Tensor]
%   t: The scaling constant. [Double]
%   g: Gamma constant in Iso Norm. [Double]
%   reg: The type of regularization to use. [Character String]
%   b: The shift in the regularization. [n x n Matrix]
%
% Outputs:
%   z: The proximal operator of tg at x. [n x n x 3 Tensor]
%      (:, :, 1) corresponds to the regularization prox.
%      (:, :, 2:3) corresponds to the iso-norm prox.

function z = prox_tg(x, t, g, reg, b)
    z = zeros(length(x(1, :, 1)), length(x(1, :, 1)), 3); % Initialize array to store output

    % Compute prox of regularization operator
    switch reg
        case 'L1' % Compute prox of L1 norm
            z(:, :, 1) = l1Prox(x(:, :, 1) - b, t) + b;
        case 'L2' % Compute prox of L2 norm
            z(:, :, 1) = l2squaredProx(x(:, :, 1) - b, t) + b; % Check this
        otherwise % Change to throw error
            error("Unrecognized regularization function specified.");
    end

    % Compute prox of Iso norm
    z(:, :, 2:3) = isoProx(x(:, :, 2:3), t*g);
end