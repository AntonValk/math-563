% primal_douglasrachford_splitting
%
% Implements non-blind image deblurring using the primal douglas rachford
% method, as described in [1]. Stops when a specified error threshold is
% met or when a maximum number of iterations have been carried out.
%
% Inputs:
%   b: The blurred image. [m x n Matrix]
%   kernel: The kernel used to blur the image. [k x k Matrix]
%   x_init: The initial guess for the deblurred image. [m x n Matrix]
%   prox_l: A method to compute the proximal operator for the regularization term. [Function Handle]
%   t: Step size. [Double]
%   g: The constant modifying the iso-norm. [Double]
%   rho: Regularization parameter. [Double]
%   k_max: Maximum number of iterations. [Integer]
%   e_t: Error threshold. [Double]
%   save: A boolean, indicating whether the image iterates should be saved. [Logical]
%   verbose: A boolean, indicating whether verbose outputs should be printed. [Logical]
%
% Outputs:
%   D: A structure containing the final image and algorithm metrics. [Struct]
%       xf: The final image. [m x n Matrix]
%       t: The time it took to run the optimization algorithm. [Double, Seconds]
%       k_end: The number of iterations ran. [Integer]
%       e_end: Error at the final iteration. [Double]
%       ek: Error at each iteration [1 x k_end Matrix]
%       xk: The image at each iteration [m x n x k_end Matrix] (Only output if save is true)
%
% Authors: Linda Hu, Cheng Shou, April Niu, Aidan Gerkis
% Date: 22-03-2024
%
% References:
%   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
%        in MATH 564 - Honours Convex Optimization.

function D = primal_douglasrachford_splitting(b, kernel, x_init, prox_l, t, g, rho, k_max, e_t, err_eval, save, verbose)
    [numRows, numCols]=size(b);
    
    % Arrays to store outputs
    errors = zeros(1, k_max);
    if save % Save images at each step only if requested
        xk = zeros(numRows, numCols, k_max);
    end

    % Initialize
    z1k = x_init;
    z2k = [mat_mult(z1k, 'K', kernel); mat_mult(z1k, 'D1', kernel); mat_mult(z1k, 'D2', kernel)];
    
    k = 1; % Current iteration
    error = e_t*10; % Current error

    tic; % Start Timer
    while error > e_t && k < k_max % Iterate until error convergence or max iterations has been exceeded
        % Split z2k into its components
        z2 = mat_split(z2k, 3);

        % Compute prox ops
        x = boxProx(z1k); % Prox of indicator
        y1 = prox_l(z2(:, :, 1), t); % Prox of regularization norm
        y2 = isoProx(z2(:, :, 2:3), t*g); % Prox of isoNorm
    
        % Compute resolvent of B
        appK = mat_mult(2*y1 - z2(:, :, 1), 'KT', kernel); % Compute the K component of the multiplication
        appD = mat_mult(2*y2 - z2(:, :, 2:3), 'DT', kernel); % Compute the D component of the multiplication (applyDTrans uses the concatenated form of the arrays)
        u = mat_mult(2*x - z1k + appK + appD, 'inv', kernel, 1); % Compute the resolvent
        
        % Compute Updates
        v = [mat_mult(u, 'K', kernel); mat_mult(u, 'D1', kernel); mat_mult(u, 'D2', kernel)];
        y = [y1; y2(:, :, 1); y2(:, :, 2)]; % Need to build y as a 2d matrix
        
        z1k = z1k + rho*(u - x);
        z2k = z2k + rho*(v - y);
        
        % Update error
        %% TODO: ERROR UPDATE <- put this in the right spot
        error = err_eval(x);

        % Save variables
        errors(k) = error;
        
        if save % Save images at each step only if requested
            xk(:, :, k) = x;
        end
        
        % Print status
        if verbose
            disp("Finished iteration " + num2str(k) + " with loss: " + error);
        end

        % Update iteration
        k = k + 1;
    end
    t_run = toc; % End Timer

    % Compile outputs
    D = struct();
    D.xf = boxProx(z1k); % Solution
    D.t = t_run; % Run time
    D.k_end = k-1; % Number of iterations
    D.e_end = 0; % Error at end
    D.ek = errors(1:D.k_end); % Error vs time
    
    if save % Save image at each iteration vs. time if requested
        D.xk = xk(:, :, 1:D.k_end);
    end
end 

