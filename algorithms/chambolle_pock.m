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
%   f: A structure containing several function handles:
%       prox_l: A method to compute the proximal operator for the regularization term. [Function Handle]
%       err_eval: A function evaluate the image error at the current iteration. [Function Handle]
%       f_loss: A function that evaluates the loss function for an image. [Function Handle]%   t: Step size. [Double]
%       early_stop: A function that evaluates the early stop condition,
%                   based on a specified tolerance. [Function Handle]
%   s: Step size. [Double]
%   g: The constant modifying the iso-norm in the problem statement. [Double]
%   k_max: Maximum number of iterations. [Integer]
%   e_t: Error threshold. [Double]
%   save: Indicates the level of saving for image iterates: [Integer]
%           0 - No image saving.
%           1 - Save images at every iterate. (WARNING: Takes lots of memory, prone to crashes)
%           2 - Sparse image saving. (Save image at every 20th Iteration)
%   ns: The step size at which to save image iterates. [Integer]%   verbose: A boolean, indicating whether verbose outputs should be printed. [Logical]
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
% Authors: Linda Hu, Cheng Shou, April Niu, Aidan Gerkis
% Date: 23-03-2024
%
% References:
%   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
%        in MATH 564 - Honours Convex Optimization.

function D = chambolle_pock(b, kernel, x_init, f, t, s, g, k_max, e_t, save, ns, verbose)
    [numRows, numCols]=size(b);
    
    % Get system information
    m = memory;
    mem_max = m.MaxPossibleArrayBytes;
    max_arr = mem_max/8; % Size of largest possible array (of double)

    % Arrays to store outputs
    errors = zeros(1, k_max);
    lossk = zeros(1, k_max);
    tk = zeros(1, k_max);
    switch save % Save images at each step only if requested
        case 0 % Do nothing, no data saved
        case 1 % Full save
            if numRows*numCols*k_max > max_arr
                error("Cannot save image iterates, not enough memory available.");
            elseif numRows*numCols*k_max*0.8 > max_arr
                warning("Image iterate array is large, will result in poor memory performance.");
            end

            xks = zeros(numRows, numCols, k_max); % Create image iterate array if possible
        case 2 % Sparse save
            if numRows*numCols*k_max/ns > max_arr
                error("Cannot save image iterates, not enough memory available.");
            elseif numRows*numCols*(k_max/ns)*0.8 > max_arr
                warning("Image iterate array is large, will result in poor memory performance.");
            end

            idx = 1; % Tracks location in the xks array
            xks = zeros(numRows, numCols, floor(k_max/ns));
        otherwise
            error("Error in function inputs: Invalid value of 'save' passed.");
    end
    
    % Initialize
    xk = x_init; % n x n matrix
    zk = x_init; % n x n matrix
    yk_1 = mat_mult(xk, 'K', kernel);
    yk_2 = mat_mult(xk, 'D', kernel);

    k = 1; % Current iteration
    error = e_t*10; % Current error
    loss_old = 0; % Initilize loss to be small
    stop = false;

    while error > e_t && k <= k_max && ~stop % Iterate until error convergence or max iterations has been exceeded
        tic; % Start Timer
        xk_old = xk; % Save previous xk
        
        % Compute Prox Op of sg*
        wk_1 = yk_1 + s*mat_mult(zk, 'K', kernel);
        wk_2 = yk_2 + s*mat_mult(zk, 'D', kernel);
        
        yk_1 = wk_1 - s*f.prox_l(wk_1/s, 1/s); % Part one of prox of sg*
        yk_2 = wk_2 - (s*g)*isoProx(wk_2/(s*g), 1/(s*g)); % Part two of prox of sg*

        % Compute Prox Op of tf
        vk = xk - t*(mat_mult(yk_1, 'KT', kernel) + mat_mult(yk_2, 'DT', kernel)); % Input to prox of f
        
        xk = boxProx(vk); % Prox of tf
       
        % Update zk
        zk = 2*xk - xk_old;
    
         % Update error & check early stop criteria
        errors(k) = f.err_eval(xk);
        lossk(k) = f.f_loss(xk);
        stop = f.early_stop(lossk(k), loss_old);
        tk(k) = toc; % End Timer

        switch save
            case 1 % Save images at each step
                xks(:, :, k) = xk;
            case 2 % Save images every ns steps
                if mod(k - 1, 20) == 0 % k - 1 ensures the first iterate is saved
                    xks(:, :, idx) = xk;
                    idx = idx + 1;
                end
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
        memory_check; % Check memory usage, should never be an issue!
    end
    t_run = sum(tk); % Time for complete deblurring process
    
    % Compile outputs
    D = struct();
    D.xf = xk; % Solution
    D.k_end = k-1; % Number of iterations
    D.t = t_run; % Run time
    D.tk = tk(1:D.k_end); % Run time for each iteration
    D.e_end = errors(D.k_end); % Error at end
    D.ek = errors(1:D.k_end); % Error vs time
    D.fk = lossk(1:D.k_end); % Loss vs iteration

    switch save % Save image at each iteration vs. time if requested
        case 1
            D.xk = xks(:, :, 1:D.k_end);
        case 2
            D.xk = xks(:, :, 1:floor(D.k_end/ns));
    end
end