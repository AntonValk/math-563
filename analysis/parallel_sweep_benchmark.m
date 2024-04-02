% parallel_sweep_benchmark.m
%
% A version of parallel_sweep that computes the simulation time and total
% speed-up, for use in benchmarking the performance of the
% parallelization. Requires that the parallel computing toolbox be
% installed.
%
% Inputs:
%   n_pool: The number of processes to allow. [Integer]
%
% Outputs:
%   sim_time: The total time to run all simulations in parallel. [Seconds]
%   speed_up: The ratio of sim time to compute time. [Unitless]
%   bytes_tf: The total number of bytes transferred. [GB]
%
% Author: Aidan Gerkis
% Date: 02-04-2024

function [sim_time, speed_up, bytes_tf] = parallel_sweep_benchmark(n_pool)
    %% Parameters
    n = 5; % Number of parameters to sweep in each direction
    k = 3; % The number of parameters to sweep, each algorithm has 3 parameters
    
    lb_1 = 1e-11; % Lower bound on parameter one
    ub_1 = 2; % Upper bound on parameter one
    
    lb_2 = 1e-11; % Lower bound on parameter two
    ub_2 = 2; % Upper bound on parameter two
    
    lb_3 = 1e-11; % Lower bound on parameter three
    ub_3 = 2; % Upper bound on parameter three
    
    % Algorithm
    alg = 'primal_dr'; % Algorithm name
    p = struct(); % Basic parameter structure to base sweep on
    p.regularization = 'L1';
    p.verbose = false;
    p.display = false;
    p.save_iters = false;
    
    %% Blur image
    I_true = imopt_scale('cameraman.jpg'); % Import image, convert to B&W with pixels in [0,1]
    
    kernel = fspecial('gaussian', [15,15], 5); % Define blur
    b = imfilter(I_true, kernel); % Apply blur
    
    noiseDensity = 0.1; % Define noise density
    b = imnoise(b, 'salt & pepper', noiseDensity); % Apply noise
    
    b = cat(3, b, I_true);
    
    %% Define params to sweep
    p1 = linspace(lb_1, ub_1, n); % Values for parameter 1
    p2 = linspace(lb_2, ub_2, n); % Values for parameter 2
    p3 = linspace(lb_3, ub_3, n); % Values for parameter 3
    
    %% Run sweep
    myPool = parpool(n_pool); % Initialize parallel environments
    idx = 1; % Index parallel calls
    
    % Initialize parallel calls
    tic;
    ticBytes(myPool); % Start counting bytes
    for i=1:n % Sweep values of parameter 1
        for j=1:n % Sweep values of parameter 2
            for k=1:n % Sweep values of parameter 3
                % Assign parameter values
                if isequal(alg, 'chambolle_pock') % Handle two cases with different parameter names
                    p.t = p1(i);
                    p.s = p2(j);
                    p.gamma = p3(k);
                else
                    p.t = p1(i);
                    p.gamma = p2(j);
                    p.rho = p3(k);
                end
    
                f_opt(idx) = parfeval(myPool, @imopt, 3, b, kernel, alg, p); % <- These start running immediately on parfeval call
                idx = idx + 1;
            end
        end
    end
    wait(f_opt);
    sim_time = toc; % Total time to run all simulations
    bytes_tf = tocBytes(myPool); % Stop counting bytes
    bytes_tf = sum(bytes_tf, 1)*1e-9; % Sum total and convert to GB

    %% Analyze Outputs
    es = zeros(n, n, n); % Store error for each parameter combination
    ts = zeros(n, n, n); % Store computation time for each parameter combination
    
    e_best = 1000; % Set initial best error to be high
    compute_time = 0; % Total time spent computing
    
    for i=1:length(f_opt) % Loop through all iterates
        % Fetch outputs, note that the first output is an index, appended by the parfeval function
        [~, ~, e_i, D_i] = fetchNext(f_opt); % <- This does not wait for all f_opt(i) to run! It processes them as they become available.
        
        if isequal(alg, 'chambolle_pock') % Extract parameters of current run
            p1_ind = find(p1 == D_i.inputs.t);
            p2_ind = find(p2 == D_i.inputs.s);
            p3_ind = find(p3 == D_i.inputs.gamma);
        else
            p1_ind = find(p1 == D_i.inputs.t);
            p2_ind = find(p2 == D_i.inputs.gamma);
            p3_ind = find(p3 == D_i.inputs.rho);
        end
        
        % Save performance metrics
        es(p1_ind, p2_ind, p3_ind) = D_i.e_end;
        ts(p1_ind, p2_ind, p3_ind) = D_i.t;
        
        compute_time = compute_time + D_i.t; % Update total compute time
    
        if e_i < e_best || i == 1 % If better performance is obtained keep this model
            D_best = D_i;
            e_best = e_i;
        end
    end
    
    speed_up = round(compute_time/sim_time, 1);
    
    %% Print best image
    imopt_display(D_best, 'Error Evolution');
    imopt_display(D_best, 'Convergence');
    if p.save_iters % If iterates were saved
        imopt_display(D_best, 'Image Iterates', 1);
    end
    
    figure('Name', 'Deblurred Image');
    imshow(D_best.xf, []);
    
    %% Delete pool
    delete(gcp('nocreate'));
end