% imopt_loss.m
% 
% Explicitly computes the loss function at a point:
%   l(x) = n_reg(Kx - b) + indicator_[0,1](x) + gamma*||Dx||_ISO
% wher n_reg is either the l1-norm or the l2-norm squared, depending on the
% problem formulation.
%
% Inputs:
%   x: The point at which to evaluate the loss function. [m x n Matrix]
%   b: The original image vector. [m x n Matrix]
%   g: Gamma, used to modify the iso-norm. [Double]
%   reg: The type of norm to use in computing the regularization term. [String]
%   kernel: The kernel of the blurring operation. [k x k Matrix]
%
% Outputs:
%   l: l(x). [Double]
%
% Author: Aidan Gerkis
% Date: 03-04-2024

function l = imopt_loss(x, b, g, reg, kernel)
    ind_s = imopt_inds(x); % Compute indicator of x on [0,1]-component wise
    n_smooth = g*imopt_iso(mat_mult(x, 'D', kernel)); % Compute gamma*||Dx||_ISO

    switch reg % Compute the regularization term
        case 'L1' % If L1 norm is used
            n_reg = imopt_l1(mat_mult(x, 'K', kernel) - b);
        case 'L2' % If L2 norm is used
            n_reg = imopt_l2(mat_mult(x, 'K', kernel) - b)^2;
        otherwise
            error('Unrecognized regularization term specified.');
    end

    l = n_reg + ind_s + n_smooth; % Compute the loss function
end