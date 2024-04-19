function p = imopt_get_default(alg, b)
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
    
    p = struct(); % Initialize
    
    % Output Parameters
    p.verbose = false; % Verbose mode disabled
    p.silent = false; % Silent mode disabled
    p.display = true; % Display plots of convergence
    p.save_iters = 0; % Disable saving of all iterates
    p.ns = 20; % Step-size at which to save iterates if sparse save is used
  
    [nr, nc] = size(b);
    p.x_init = zeros(nr, nc); % Initial Starting Point - 0 by default
    p.x_true = NaN(nr, nc); % The true image, all NaNs indicates it was not passed

    % Add global algorithm parameters
    p.max_iter = 500; % Maximum number of iterations
    p.e_t = 0.1; % Error threshold
    p.tol = 0; % Loss threshold
    p.metric = 'rmse'; % Error metric to use when evaluating convergence
    p.regularization = 'L1'; % Regularization type
    
    % Add algorithm specific parameters
    switch alg 
        case 'primal_dr'
            p.t = 0.1; % Step size
            p.gamma = 0.01; % Gamma
            p.rho = 1.05; % Regularization constant
        case 'primaldual_dr'
            p.t = 1.75; % Step size
            p.gamma = 0.05; % Gamma
            p.rho = 1.75; % Regularization constant
        case 'admm'
            p.t = 0.7407; % Step size
            p.gamma = 1e-15; % Gamma
            p.rho = 1.4815;%1.34815; % Regularization constant
        case 'chambolle_pock'
            p.t = 0.1; % Step size
            p.s = 0.4072; % Step size
            p.gamma = 1e-15; % Gamma
        otherwise
            error("Unrecognized algorithm specified.");
    end

    % Parameters for saving image
    p.im_name = ""; % Name of the image
    p.dir = convertCharsToStrings(pwd); % Default path is current folder
end