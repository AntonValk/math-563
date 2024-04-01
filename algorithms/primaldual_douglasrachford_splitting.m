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
%   prox_g: The proximal operator computation for g(x), the regularization
%           terms. [Function Handle]
%   t: Step size. [Double]
%   g: The constant modifying the iso-norm in the problem statement. [Double]
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
% Date: 23-03-2024
%
% References:
%   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
%        in MATH 564 - Honours Convex Optimization.

function  D = primaldual_douglasrachford_splitting(b, kernel, x_init, prox_g, t, g, rho, k_max, e_t, save, verbose)
    [numRows, numCols]=size(b);
    
    % Arrays to store outputs
    errors = zeros(1, k_max);
    if save % Save images at each step only if requested
        xks = zeros(numRows, numCols, k_max);
    end
    
    % Initialize
    pk = x_init; % n x n matrix
    qk = [mat_mult(pk, 'K', kernel); mat_mult(pk, 'D1', kernel); mat_mult(pk, 'D2', kernel)]; % 3n x n matrix
    
    k = 1; % Current iteration
    error = e_t*10; % Current error
    
    tic; % Start Timer
    while error > e_t && k < k_max % Iterate until error convergence or max iterations has been exceeded
        q = mat_split(qk, 3); % Convert qk to a 3D tensor
        q_prox = q;
        q_prox(:, :, 1) = q_prox(:, :, 1)/t;
        q_prox(:, :, 2:3) = q_prox(:, :, 2:3)/(t*g);
    
        % Compute Prox Ops
        xk = boxProx(pk);   %x is n by n matrix
    
        zk_conj = prox_g(q_prox, 1/t); % Prox op of g
    
        z1 = q(:, :, 1) - t*zk_conj(:, :, 1); % Prox op of g* (regularization part)
        z2 = q(:, :, 2:3) - (t*g)*zk_conj(:, :, 2:3); % Prox op of g* (iso-norm part)
    
        %z1 = q(:, :, 1) - t*l1Prox(q(:, :, 1)/t - b, 1/t) - b; % Equation 8 in reference
        %z2 = q(:, :, 2:3) - (t*g)*isoProx(q(:, :, 2:3)/(t*g), 1/(t*g)); % Why using g
        z21 = z2(:,:,1);
        z22 = z2(:,:,2);
        zk = [z1; z2(:, :, 1); z2(:, :, 2)];
    
        % Compute Resolvent of B (pg 7 in reference) <- TODO: Double check
        vec0  =[2*xk - pk; 2*z1 - q(:,:,1); 2*z21-q(:,:,2);2*z22-q(:,:,3)];
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
    D.xf = boxProx(pk); % Solution
    D.t = t_run; % Run time
    D.k_end = k-1; % Number of iterations
    D.e_end = 0; % Error at end
    D.ek = errors(1:D.k_end); % Error vs time
    
    if save % Save image at each iteration vs. time if requested
        D.xk = xks(:, :, 1:D.k_end);
    end
end