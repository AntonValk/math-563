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
%   prox_g: The proximal operator computation for g(x), the regularization
%           terms. [Function Handle]
%   t: Step size. [Double]
%   s: Step size. [Double]
%   g: The constant modifying the iso-norm in the problem statement. [Double]
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
% Date: 23-03-2024
%
% References:
%   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
%        in MATH 564 - Honours Convex Optimization.

function D = chambolle_pock(b, kernel, x_init, prox_g, t, s, g, k_max, e_t, save, verbose)
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
        
        wk_prox = wk;
        wk_prox(:, :, 1) = wk_prox(:, :, 1)/s; % Scale wk for input to prox (moreau decomposition)
        wk_prox(:, :, 2:3) = wk_prox(:, :, 2:3)/(s*g); % Scale wk for input to prox (moreau decomposition)
    
        % Prox of sg*
        yk_conj = prox_g(wk_prox, 1/s); % Prox op of g
        
        yk1 = wk(:, :, 1) - s*yk_conj(:, :, 1); % Prox op of sg* (regularization part)
        yk2 = wk(:, :, 2:3) - (s*g)*yk_conj(:, :, 2:3); % Prox op of g* (iso-norm part)
    
        %yk1 = wk(:, :, 1) - s*l1Prox(wk(:, :, 1)/s - b, 1/s) - b; % Part one of prox of sg*
        %yk2 = wk(:, :, 2:3) - s*isoProx(wk(:, :, 2:3)/s, g/s); % Part two of prox of sg*
    
        yk = [yk1; yk2(:, :, 1); yk2(:, :, 2)]; % Compile components of prox of sg*
    
        % Prox of tf
        xk = boxProx(vk);
    
        % Update zk
        zk = 2*xk - xk_old;
    
        % Update error
        %% TODO: ERROR UPDATE
    
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
    D.e_end = 0; % Error at end
    D.ek = errors(1:D.k_end); % Error vs time
    
    if save % Save image at each iteration vs. time if requested
        D.xk = xks(:, :, 1:D.k_end);
    end
end