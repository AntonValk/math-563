function D = primal_douglasrachford_splitting(b, kernel, x_init, f, t, g, rho, k_max, e_t, save, ns, verbose)
    % primal_douglasrachford_splitting.m
    %
    % Implements non-blind image deblurring using the primal douglas rachford
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
    %   save: Indicates the level of saving for image iterates: [Integer]
    %           0 - No image saving.
    %           1 - Save images at every iterate. (WARNING: Takes lots of memory, prone to crashes)
    %           2 - Sparse image saving. (Save image at every 20th Iteration)
    %   ns: The step size at which to save image iterates. [Integer]
    %   verbose: A boolean, indicating whether verbose outputs should be printed. [Logical]
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
    % Date: 22-03-2024
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
    z1k = x_init;
    z2k_1 = mat_mult(z1k, 'K', kernel); % Store z2 as two parts: m x n component
    z2k_2 = mat_mult(z1k, 'D', kernel); % Store z2 as two parts: m x n x 2 component

    k = 1; % Current iteration
    err = e_t*10; % Current error
    loss_old = 0; % Initilize loss to be small
    stop = false;

    while err > e_t && k <= k_max && ~stop % Iterate until error convergence or max iterations has been exceeded
        tic; % Start Timer

        % Compute prox ops
        x = boxProx(z1k); % Prox of indicator
        y1 = f.prox_l(z2k_1, t); % Prox of regularization norm
        y2 = isoProx(z2k_2, t*g); % Prox of isoNorm
    
        % Compute resolvent of B
        appK = mat_mult(2*y1 - z2k_1, 'KT', kernel); % Compute the K component of the multiplication
        appD = mat_mult(2*y2 - z2k_2, 'DT', kernel); % Compute the D component of the multiplication (applyDTrans uses the concatenated form of the arrays)
        u = mat_mult(2*x - z1k + appK + appD, 'inv', kernel, 1); % Compute the resolvent

        v_1 = mat_mult(u, 'K', kernel); % Compute Au = [Ku; Du]
        v_2 = mat_mult(u, 'D', kernel);
        
        % Compute Updates
        z1k = z1k + rho*(u - x);
        z2k_1 = z2k_1 + rho*(v_1 - y1);
        z2k_2 =  z2k_2 + rho*(v_2 - y2);
        
        % Update error & check early stop criteria
        errors(k) = f.err_eval(x);
        lossk(k) = f.f_loss(x);
        stop = f.early_stop(lossk(k), loss_old);
        tk(k) = toc; % End Timer
        
        switch save
            case 1 % Save images at each step
                xks(:, :, k) = x;
            case 2 % Save images every ns steps
                if mod(k - 1, 20) == 0 % k - 1 ensures the first iterate is saved
                    xks(:, :, idx) = x;
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
    D.xf = x; % Solution
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

    if D.k_end < k_max && verbose % Communicate early stop
        fprintf("Early stop at iteration %d.\n", D.k_end);
    end
end