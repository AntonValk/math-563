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
%   f: A structure containing several function handles:
%       prox_l: A method to compute the proximal operator for the regularization term. [Function Handle]
%       err_eval: A function evaluate the image error at the current iteration. [Function Handle]
%       f_loss: A function that evaluates the loss function for an image. [Function Handle]
%       early_stop: A function that evaluates the early stop condition,
%                   based on a specified tolerance. [Function Handle]
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
%       tk: The time it took to compute iteration k. [Double, Seconds]
%       k_end: The number of iterations ran. [Integer]
%       e_end: Error at the final iteration. [Double]
%       ek: Error at each iteration [1 x k_end Matrix]
%       fk: Loss at each iteration [1 x k_end Matrix]
%       xk: The image at each iteration [m x n x k_end Matrix] (Only output if save is true)
%
% Authors: Linda Hu, Cheng Shou, April Niu
% Date: 22-03-2024
%
% References:
%   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
%        in MATH 564 - Honours Convex Optimization.

function D = admm(b, kernel, x_init, f, t, rho, g, k_max, e_t, save, verbose)
    [numRows, numCols]=size(b);
    
    % Arrays to store outputs
    errors = zeros(1, k_max);
    lossk = zeros(1, k_max);
    tk = zeros(1, k_max);
    if save % Save images at each step only if requested
        xks = zeros(numRows, numCols, k_max);
    end
    
    % Initialize
    xk = x_init; % m x n matrix
    uk = x_init; % m x n matrix
    wk = x_init; % m x n matrix
    y1 = mat_mult(xk, 'K', kernel); % m x n matrix
    y2 = mat_mult(xk, 'D', kernel); % m x n x 2 matrix
    z1 = mat_mult(xk, 'K', kernel); % m x n matrix
    z2 = mat_mult(xk, 'D', kernel); % m x n x 2 matrix
    
    k = 1; % Current iteration
    error = e_t*10; % Current error
    loss_old = 0; % Initilize loss to be small
    stop = false;

    while error > e_t && k <= k_max && ~stop % Iterate until error convergence or max iterations has been exceeded
        tic; % Start Timer
        % Compute inputs to prox op calculations
        atz = mat_mult(z1, 'KT', kernel) + mat_mult(z2, 'DT', kernel); %(A^t z)
        aty = mat_mult(y1, 'KT', kernel) + mat_mult(y2, 'DT', kernel); %(A^t y)
        xk = mat_mult(uk + aty - (1/t)*(wk + atz), 'inv', kernel, 1); % Update image

        ax_1 = mat_mult(xk, 'K', kernel); % Compute Ax - m x n
        ax_2 = mat_mult(xk, 'D', kernel); % Compute Ax - m x n x 2

        % Compute prox ops
        uk = boxProx(rho*xk + (1-rho)*uk + wk/t); % Prox of f/t
        
        y1 = f.prox_l(rho*ax_1 + (1-rho)*y1 + z1/t, 1/t); % Part 1 of prox g/t
        y_aux = isoProx(rho*ax_2 + (1-rho)*y2 + z2/t, g/t); % Part 2 of prox g/t
    
        % Compute Updates
        wk = wk + t*(xk-uk);
        z1 = z1 + t * (mat_mult(xk, 'K', kernel) - y1); % Component-wise update of zk - m x n
        z2 = z2 + t * (mat_mult(xk, 'D', kernel) - y_aux); % Component-wise update of zk - m x n x 2
        
        % Update error & check early stop criteria
        errors(k) = f.err_eval(uk);
        lossk(k) = f.f_loss(uk);
        stop = f.early_stop(lossk(k), loss_old);
        tk(k) = toc; % End Timer

        if save % Save images at each step only if requested
            xks(:, :, k) = uk;
        end
    
        % Print status
        if verbose
            disp("Finished iteration " + num2str(k) + " in " + num2str(tk(k)) + "s with: ");
            disp("  error: " + num2str(errors(k)));
            disp("  loss: " + num2str(lossk(k)));
        end

        % Update iteration
        loss_old = lossk(k);
        k = k + 1;
    end
    
    t_run = sum(tk); % Time for complete deblurring process
    
    % Compile outputs
    D = struct();
    D.xf = boxProx(xk); % Solution. Perform additional boxProx to enforce constraints
    D.k_end = k-1; % Number of iterations
    D.t = t_run; % Run time
    D.tk = tk(1:D.k_end); % Run time for each iteration
    D.e_end = errors(D.k_end); % Error at end
    D.ek = errors(1:D.k_end); % Error vs iteration
    D.fk = lossk(1:D.k_end); % Loss vs iteration

    if save % Save image at each iteration vs. time if requested
        D.xk = xks(:, :, 1:D.k_end);
    end
end