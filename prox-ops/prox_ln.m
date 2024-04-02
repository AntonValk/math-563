% prox_ln.m
%
% Evaluates the prox operator of the L1 or L2 norm, as specified by the
% user.
%
% Inputs:
%   x: The point at which to evaluate the prox. [m x n Matrix]
%   t: The scaling constant. [Double]
%   norm: The type of norm to use. Can be one of: [Character String]
%           - 'L1'
%           - 'L2'
%   b: The shift in the regularization. [m x n Matrix]
%
% Outputs:
%   z: The proximal operator of the specified norm at x. [m x n Matrix]

function z = prox_ln(x, t, norm, b)
    switch norm % Compute prox of norm
        case 'L1' % Compute prox of L1 norm
            z = l1Prox(x(:, :, 1) - b, t) + b;
        case 'L2' % Compute prox of L2 norm
            z = l2squaredProx(x(:, :, 1) - b, t) + b; % Check this
        otherwise % Change to throw error
            disp("Error in prox_ln: Unrecognized norm specified.");
            return
    end
end