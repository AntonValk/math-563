% primaldual_douglasrachford_splitting
%
% Implements non-blind image deblurring using the primal-dual douglas rachford
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
%       fk: Loss at each iteration [1 x k_end Matrix]
%       xk: The image at each iteration [m x n x k_end Matrix] (Only output if save is true)
%
% Authors: Linda Hu, Cheng Shou, April Niu, Aidan Gerkis
% Date: 23-03-2024
%
% References:
%   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
%        in MATH 564 - Honours Convex Optimization.

function  D = primaldual_douglasrachford_splitting(b, kernel, x_init, f, t, g, rho, k_max, e_t, save, verbose)
    [numRows, numCols]=size(b);
    
    % Arrays to store outputs
    errors = zeros(1, k_max);
    lossk = zeros(1, k_max);
    if save % Save images at each step only if requested
        xks = zeros(numRows, numCols, k_max);
    end
    
    % Initialize
    pk = x_init; % n x n matrix
    qk = [mat_mult(pk, 'K', kernel); mat_mult(pk, 'D1', kernel); mat_mult(pk, 'D2', kernel)]; % 3n x n matrix
    
    k = 1; % Current iteration
    error = e_t*10; % Current error
    loss_old = 0; % Initilize loss to be small
    stop = false;

    tic; % Start Timer
    while error > e_t && k < k_max && ~stop % Iterate until error convergence or max iterations has been exceeded
        q = mat_split(qk, 3); % Convert qk to a 3D tensor
    
        % Compute Prox Ops
        xk = boxProx(pk);   %x is n by n matrix
        
        z1 = q(:, :, 1) - t*f.prox_l(q(:, :, 1)/t, 1/t); % Part 1 of prox_sg*
        z2 = q(:, :, 2:3) - (t*g)*isoProx(q(:, :, 2:3)/(t*g), 1/(t*g)); % Part 2 of prox_sg*
        zk = [z1; z2(:, :, 1); z2(:, :, 2)]; % Compile prox_sg*
    
        % Compute Resolvent of B (pg 7 in reference)
        vec0  =[2*xk - pk; 2*z1 - q(:,:,1); 2*z2(:,:,1) - q(:,:,2); 2*z2(:,:,2) - q(:,:,3)];
        vec = mat_split(vec0, 4); % Extract matrices corresponding to n x n blocks
    
        a = (vec(:, :, 1)) - t*mat_mult(vec(:, :, 2), 'KT', kernel) - t*mat_mult(vec(:, :, 3), 'D1T', kernel) - t*mat_mult(vec(:, :, 4), 'D2T', kernel); % [I, -tA']*vec
        b = mat_mult(a, 'inv', kernel, t); % (I + t^2A^TA)^-1 * [I, -tA']*vec
        c = [eye(numRows)*b; t*mat_mult(b, 'K', kernel); t*mat_mult(b, 'D1', kernel); t*mat_mult(b, 'D2', kernel)]; % [I; tA]*(I + t^2A^TA)^-1 * [I, -tA']*vec
    
        res = [zeros(numRows,numCols); vec(:, :, 2); vec(:, :, 3); vec(:, :, 4)] + c;
    
        wk = res(1:numRows, :); % In n x n
        vk = res((numRows + 1):end, :); % In 3n x n
    
        % Perform update
        pk = pk + rho*(wk - xk); % n x n matrix
        qk = qk + rho*(vk - zk); % 3n x n matrix
        
        % Update error & check early stop criteria
        errors(k) = f.err_eval(xk);
        lossk(k) = f.f_loss(xk);
        stop = f.early_stop(lossk(k), loss_old);

        if save % Save images at each step only if requested
            xks(:, :, k) = xk;
        end
    
        % Print status
        if verbose
            disp("Finished iteration " + num2str(k) + " with: ");
            disp("  error: " + num2str(errors(k)));
            disp("  loss: " + num2str(lossk(k)));
        end

        % Update iteration
        loss_old = lossk(k);
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
    D.fk = lossk(1:D.k_end); % Loss vs iteration

    if save % Save image at each iteration vs. time if requested
        D.xk = xks(:, :, 1:D.k_end);
    end
end