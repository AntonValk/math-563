
function D = admm(b, kernel, x_init, prox_l, t, rho, g, k_max, e_t, save, verbose)
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
        atz = mat_mult(z1, 'KT', kernel) + mat_mult(z2, 'D1T', kernel) + mat_mult(z3, 'D2T', kernel); %(A^t z)
        aty = mat_mult(y1, 'KT', kernel) + mat_mult(y2, 'D1T', kernel) + mat_mult(y3, 'D2T', kernel); %(A^t y)
        xk = mat_mult(uk + aty - (1/t)*(wk + atz), 'inv', kernel, 1); %
    
        % Compute prox ops
        uk = boxProx(rho*xk + (1-rho)*uk + wk/t);
        w = mat_split([mat_mult(xk, 'K', kernel); mat_mult(xk, 'D1', kernel); mat_mult(xk, 'D2', kernel)], 3);
        yk = mat_split([y1; y2; y3], 3);
        zk = mat_split([z1; z2; z3], 3);
        % unsure about the +/-b and */t for the two prox operators.
        y1 = prox_l(rho* w(:,:,1) + (1-rho)*y1+z1/t,1/t); % From Equation 8 in the reference
        y_aux = isoProx(rho*w(:,:,2:3) + (1-rho)*yk(:,:,2:3) + zk(:,:,2:3)/t, g/t);
        %yk = [y1; y_aux(:,:,1); y_aux(:,:,2)];
    
        % Compute Updates
        wk = wk + t*(xk-uk);
        z1 = z1 + t * (mat_mult(xk, 'K', kernel) - y1);
        z2 = z2 + t * (mat_mult(xk, 'D1', kernel) - y_aux(:,:,1));
        z3 = z3 + t * (mat_mult(xk, 'D2', kernel) - y_aux(:,:,2));
    
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
    atz = mat_mult(z1, 'KT', kernel) + mat_mult(z2, 'D1T', kernel) + mat_mult(z3, 'D2T', kernel); %(A^t z)
    aty = mat_mult(y1, 'KT', kernel) + mat_mult(y2, 'D1T', kernel) + mat_mult(y3, 'D2T', kernel); %(A^t y)
    xSol= mat_mult(uk + aty - (1/t)*(wk + atz), 'inv', kernel, 1);
    
    t_run = toc; % End Timer
    
    % Compile outputs
    D = struct();
    D.xf = xSol; % Solution
    D.t = t_run; % Run time
    D.k_end = k-1; % Number of iterations
    D.e_end = 0; % Error at end
    D.ek = errors(1:D.k_end); % Error vs time
    
    if save % Save image at each iteration vs. time if requested
        D.xk = xks(:, :, 1:D.k_end);
    end
end