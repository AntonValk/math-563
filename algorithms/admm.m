% admm
%
% Implements non-blind image deblurring using the admm
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
% Authors: Linda Hu, Cheng Shou, April Niu
% Date: 22-03-2024
%
% References:
%   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
%        in MATH 564 - Honours Convex Optimization.

function D = admm(b, kernel, x_init, prox_l, t, rho, g, k_max, e_t, err_eval, save, verbose)
    [numRows, numCols]=size(b);
    
    % Arrays to store outputs
    errors = zeros(1, k_max);
    if save % Save images at each step only if requested
        xks = zeros(numRows, numCols, k_max);
    end
    
    % Initialize
    xk = x_init; % m x n matrix
    uk = x_init; % m x n matrix
    wk = x_init; % m x n matrix
    y1 = mat_mult(xk, 'K', kernel);  % 3m x n matrix
    y2 = mat_mult(xk, 'D1', kernel); % 3m x n matrix
    y3 = mat_mult(xk, 'D2', kernel); % 3m x n matrix
    z1 = mat_mult(xk, 'K', kernel); % 3m x n matrix
    z2 = mat_mult(xk, 'D1', kernel); % 3m x n matrix
    z3 = mat_mult(xk, 'D2', kernel); % 3m x n matrix
    
    k = 1; % Current iteration
    error = e_t*10; % Current error
    
    tic; % Start Timer
    while error > e_t && k < k_max % Iterate until error convergence or max iterations has been exceeded
        % Compute inputs to prox op calculations
        atz = mat_mult(z1, 'KT', kernel) + mat_mult(z2, 'D1T', kernel) + mat_mult(z3, 'D2T', kernel); %(A^t z)
        aty = mat_mult(y1, 'KT', kernel) + mat_mult(y2, 'D1T', kernel) + mat_mult(y3, 'D2T', kernel); %(A^t y)
        xk = mat_mult(uk + aty - (1/t)*(wk + atz), 'inv', kernel, 1); % 

        ax = mat_split([mat_mult(xk, 'K', kernel); mat_mult(xk, 'D1', kernel); mat_mult(xk, 'D2', kernel)], 3); % Compute Ax and convert 2D -> 3D
        yk = mat_split([y1; y2; y3], 3); % Convert 2D -> 3D
        zk = mat_split([z1; z2; z3], 3); % Convert 2D -> 3D

        % Compute prox ops
        uk = boxProx(rho*xk + (1-rho)*uk + wk/t); % Prox of f/t
        
        y1 = prox_l(rho*ax(:,:,1) + (1-rho)*y1 + z1/t, 1/t); % Part 1 of prox g/t
        y_aux = isoProx(rho*ax(:,:,2:3) + (1-rho)*yk(:,:,2:3) + zk(:,:,2:3)/t, g/t); % Part 2 of prox g/t
    
        % Compute Updates
        wk = wk + t*(xk-uk);
        z1 = z1 + t * (mat_mult(xk, 'K', kernel) - y1); % Component-wise update of zk
        z2 = z2 + t * (mat_mult(xk, 'D1', kernel) - y_aux(:,:,1));
        z3 = z3 + t * (mat_mult(xk, 'D2', kernel) - y_aux(:,:,2));
    
        % Update error
        error = err_eval(xk);
    
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