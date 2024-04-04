% get_default.m
%
% Returns the default parameter settings for the specified algorithm.
%
% Inputs:
%   alg: The algorithm for which to return the default parameters. [String]
%
% Outputs:
%   p: The parameter structure
%
% Author: Aidan Gerkis
% Date: 31-3-2024

function p = get_default(alg, b)
    p = struct(); % Initialize
    
    % Output Parameters
    p.verbose = 0; % Verbose mode disabled
    p.display = 1; % Display plots of convergence
    p.save_iters = 0; % Disable saving of all iterates
    
    % Initial Starting Point - 0 by default
    [nr, nc] = size(b);
    p.x_init = zeros(nr, nc);

    % Add global algorithm parameters
    p.max_iter = 500; % Maximum number of iterations
    p.e_t = 0.1; % Error threshold
    p.e_meas = 'delta'; % Error metric to use when evaluating convergence
    p.regularization = 'L1'; % Regularization type
    p.metric = 'rmse'; % Error metric to use
    
    % Add algorithm specific parameters
    switch alg 
        case 'primal_dr'
            p.t = 0.25; % Step size
            p.gamma = 0.15; % Gamma
            p.rho = 0.85; % Regularization constant
        case 'primaldual_dr'
            p.t = 0.5; % Step size
            p.gamma = 0.05; % Gamma
            p.rho = 1.5; % Regularization constant
        case 'admm'
            p.t = 1; % Step size
            p.gamma = 0.1; % Gamma
            p.rho = 0.5; % Regularization constant
        case 'chambolle_pock'
            p.t = 0.35; % Step size
            p.s = 0.35; % Step size
            p.gamma = 0.25; % Gamma
        otherwise
            error("Unrecognized algorithm specified.");
    end
end