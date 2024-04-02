% chambolle_pock.m
%
% Implements non-blind image deblurring using the chambolle-pock
% method, as described in [1]. Stops when a specified error threshold is
% met or when a maximum number of iterations have been carried out.
%
% Inputs:
%   b: The blurred image. [m x n Matrix]
%   kernel: The kernel used to blur the image. [k x k Matrix]
%   x_init: The initial guess for the deblurred image. [m x n Matrix]
%   prox_l: A method to compute the proximal operator for the regularization term. [Function Handle]
%   t: Step size. [Double]
%   s: Step size. [Double]
%   g: The constant modifying the iso-norm in the problem statement. [Double]
%   k_max: Maximum number of iterations. [Integer]
%   e_t: Error threshold. [Double]
%   err_eval: A function evaluate the image error at the current iteration. [Function Handle]
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
% Date: 23-03-2024
%
% References:
%   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
%        in MATH 564 - Honours Convex Optimization.

function D = chambolle_pock(b, kernel, x_init, prox_l, t, s, g, k_max, e_t, err_eval, save, verbose)
    [numRows, numCols]=size(b);
    
    % Arrays to store outputs
    errors = zeros(1, k_max);
    if save % Save images at each step only if requested
        xks = zeros(numRows, numCols, k_max);
    end
    
    
    % Initialize
    xk = x_init; % n x n matrix
    zk = x_init; % n x n matrix
    yk = [mat_mult(xk, 'K', kernel); mat_mult(xk, 'D1', kernel); mat_mult(xk, 'D2', kernel)]; % 3n x n matrix
    
    k = 1; % Current iteration
    error = e_t*10; % Current error
    
    tic; % Start Timer
    while error > e_t && k < k_max % Iterate until error convergence or max iterations has been exceeded
        xk_old = xk; % Save previous xk
    
        % Compute Prox Ops
        wk = yk + s*[mat_mult(zk, 'K', kernel); mat_mult(zk, 'D1', kernel); mat_mult(zk, 'D2', kernel)]; % Input to prox of g
        wk = mat_split(wk, 3); % Convert from 2D -> 3D tensor
    
        yk = mat_split(yk, 3); % Convert from 2D -> 3D tensor
        vk = xk - t*(mat_mult(yk(:, :, 1), 'KT', kernel) + mat_mult(yk(:, :, 2), 'D1T', kernel) + mat_mult(yk(:, :, 3), 'D2T', kernel)); % Input to prox of f
        
        xk = boxProx(vk); % Prox of tf
        yk1 = wk(:, :, 1) - s*prox_l(wk(:, :, 1)/s, 1/s); % Part one of prox of sg*
        yk2 = wk(:, :, 2:3) - (s*g)*isoProx(wk(:, :, 2:3)/(s*g), 1/(s*g)); % Part two of prox of sg*
    
        yk = [yk1; yk2(:, :, 1); yk2(:, :, 2)]; % Compile components of prox of sg*
    
        % Update zk
        zk = 2*xk - xk_old;
    
        % Update error
        error = err_eval(x);
    
        % Save variables
        errors(k) = error;
    
        if save % Save images at each step only if requested
            xks(:, :, k) = xk;
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
    D.xf = xk; % Solution
    D.t = t_run; % Run time
    D.k_end = k-1; % Number of iterations
    D.e_end = errors(k-1); % Error at end
    D.ek = errors(1:D.k_end); % Error vs time
    
    if save % Save image at each iteration vs. time if requested
        D.xk = xks(:, :, 1:D.k_end);
    end
end