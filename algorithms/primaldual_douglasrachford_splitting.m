function  D = primaldual_douglasrachford_splitting(b, kernel, x_init, f, t, g, rho, k_max, e_t, save, ns, verbose)
    % primaldual_douglasrachford_splitting.m
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
    %       e_end: Error at the final iteration. Only output if the true image is passed. [Double]
    %       ek: Error at each iteration. Only output if the true image is passed. [1 x k_end Matrix]
    %       fk: Loss at each iteration [1 x k_end Matrix]
    %       xk: The image at each iteration [m x n x k_end Matrix] (Only output if save is true)
    %
    % Authors: Linda Hu, Cheng Shou, April Niu, Aidan Gerkis
    % Date: 23-03-2024
    %
    % References:
    %   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
    %        in MATH 564 - Honours Convex Optimization.
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
    pk = x_init; % n x n matrix
    qk_1 = mat_mult(pk, 'K', kernel); % Store qk as two components: m x n component
    qk_2 = cat(3, mat_mult(pk, 'D1', kernel), mat_mult(pk, 'D2', kernel)); % Store qk as two components: m x n x 2 component

    k = 1; % Current iteration
    err = e_t*10; % Current error
    loss_old = 0; % Initilize loss to be small
    stop = false;

    while err > e_t && k <= k_max && ~stop % Iterate until error convergence or max iterations has been exceeded
        tic; % Start Timer
    
        % Compute Prox Ops
        xk = boxProx(pk); % Prox_tf
        
        z1 = qk_1 - t*f.prox_l(qk_1/t, 1/t); % Part 1 of prox_sg*
        z2 = qk_2 - (t*g)*isoProx(qk_2/(t*g), 1/(t*g)); % Part 2 of prox_sg*
    
        % Compute Resolvent of B (pg 7 in reference)
        vec_1 = 2*xk - pk; % m x n
        vec_2 = 2*z1 - qk_1; % m x n
        vec_3 = 2*z2 - qk_2; % m x n x 2

        a = (vec_1) - t*mat_mult(vec_2, 'KT', kernel) - t*mat_mult(vec_3, 'DT', kernel); % [I, -tA']*vec
        b = mat_mult(a, 'inv', kernel, t); % (I + t^2A^TA)^-1 * [I, -tA']*vec
        wk = eye(numRows)*b; % [0, 0; 0, I]*vec + [I; tA]*(I + t^2A^TA)^-1 * [I, -tA']*vec [m x n]
        vk_1 = vec_2 + t*mat_mult(b, 'K', kernel); % [0, 0; 0, I]*vec + [I; tA]*(I + t^2A^TA)^-1 * [I, -tA']*vec [m x n]
        vk_2 = vec_3 + t*mat_mult(b, 'D', kernel); % [0, 0; 0, I]*vec + [I; tA]*(I + t^2A^TA)^-1 * [I, -tA']*vec [m x n x 2]
    
        % Perform update
        pk = pk + rho*(wk - xk); % n x n matrix
        qk_1 = qk_1 + rho*(vk_1 - z1);
        qk_2 = qk_2 + rho*(vk_2 - z2);
        
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
            if ~isnan(errors(k)) % Display error if ground truth was passed
                disp("  error: " + num2str(errors(k)));
            end
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
    if ~anynan(errors) % If true image was passed
        D.e_end = errors(D.k_end); % Error at end
        D.ek = errors(1:D.k_end); % Error vs time
    end
    D.fk = lossk(1:D.k_end); % Loss vs iteration

    switch save % Save image at each iteration vs. time if requested
        case 1
            D.xk = xks(:, :, 1:D.k_end);
        case 2
            D.xk = xks(:, :, 1:floor(D.k_end/ns));
    end
end