function D = primal_douglasrachford_splitting(b, kernel, x_init, prox_g, t, rho, k_max, e_t, save)
    [numRows, numCols]=size(b);
    
    % Arrays to store outputs
    errors = zeros(1, k_max);
    if save % Save images at each step only if requested
        xk = zeros(numRows, numCols, k_max);
    end

    % Initialize
    z1k = x_init;
    z2k = [mat_mult(z1k, 'K', kernel); mat_mult(z1k, 'D1', kernel); mat_mult(z1k, 'D2', kernel)];
    
    k = 1;
    error = e_t*10;

    tic; % Start Timer
    while error > e_t && k < k_max
        % Split z2k into its components
        z2 = mat_split(z2k, 3);
        %z21 = z2k(1:numRows, :);
        %z22 = cat(3, z2k((numRows + 1):2*numRows, :), z2k((2*numRows + 1):end, :)); % Change z2 from 2D (2n*n)->3D (n*n*2)
        
        % Compute prox ops
        x = boxProx(z1k);
        y = prox_g(z2, t);
        y1 = y(:, :, 1);
        y2 = y(:, :, 2:3);
    
        %y1 = l1Prox(z2(:, :, 1) - b, t) + b; % From Equation 8 in the reference
        %y2 = isoProx(z2(:, :, 2:3), t*g);
    
        % Compute resolvent of B
        appK = mat_mult(2*y1 - z2(:, :, 1), 'KT', kernel); % Compute the K component of the multiplication
        appD = mat_mult(2*y2 - z2(:, :, 2:3), 'DT', kernel); % Compute the D component of the multiplication (applyDTrans uses the concatenated form of the arrays)
        u = mat_mult(2*x - z1k + appK + appD, 'inv', kernel, 1); % Compute the resolvent
        
        % Compute Updates
        v = [mat_mult(u, 'K', kernel); mat_mult(u, 'D1', kernel); mat_mult(u, 'D2', kernel)];
        y = [y1; y2(:, :, 1); y2(:, :, 2)]; % Need to build y as a 2d matrix
        
        z1k = z1k + rho*(u - x);
        z2k = z2k + rho*(v - y);
        
        % Save variables
        errors(k) = error;
        
        if save % Save images at each step only if requested
            xk(:, :, k) = x;
        end

        % Update iteration
        k = k + 1;
    end
    t_run = toc; % End Timer

    % Compile outputs
    D = struct();
    D.xf = boxProx(z1k); % Solution
    D.t = t_run; % Run time
    D.k_end = k-1; % Number of iterations
    D.e_end = 0; % Error at end
    D.ek = errors(1:D.k_end); % Error vs time
    
    if save % Save image at each iteration vs. time if requested
        D.xk = xk(:, :, 1:D.k_end);
    end
end 

