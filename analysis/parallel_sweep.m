% parallel_sweep.m
%
% Performs a sweep of the three model hyperparameters, utilizing MATLAB's
% parallelization toolbox. Requires that the parallel computing toolbox be
% installed.
%
% Author: Aidan Gerkis
% Date:02-04-2024

clear; clc; close all;

%% Parameters
n = 25; % Number of parameters to sweep in each direction
k = 3; % The number of parameters to sweep, each algorithm has 3 parameters

lb_1 = 1e-11; % Lower bound on parameter one
ub_1 = 2; % Upper bound on parameter one

lb_2 = 1e-11; % Lower bound on parameter two
ub_2 = 2; % Upper bound on parameter two

lb_3 = 1e-11; % Lower bound on parameter three
ub_3 = 2; % Upper bound on parameter three

fn_out = "n25-largesweep";

% Parallelization
n_pool = 16; % Number of processes to allow. Chosen from benchmarking on a 24 core computer
%n_pool = 8; % Number of processes to allow. Use this for your PCs! Or maybe even lower, check your system info.

% Algorithm
alg = 'primaldual_dr'; % Algorithm name
p = struct(); % Basic parameter structure to base sweep on
p.regularization = 'L1';
p.verbose = false;
p.display = false;
p.save_iters = false;

%% Blur image
I_true = imopt_scale('cameraman.jpg'); % Import image, convert to B&W with pixels in [0,1]

kernel = fspecial('gaussian', [9,9], 4); % Define blur
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

%% Analyze Outputs
fs = zeros(n, n, n); % Store loss for each parameter combination
es = zeros(n, n, n); % Store error for each parameter combination
ts = zeros(n, n, n); % Store computation time for each parameter combination

f_best = Inf; % Set initial best loss to be high

for i=1:length(f_opt) % Loop through all iterates
    % Fetch outputs, note that the first output is an index, appended by the parfeval function
    [~, x, e_i, D_i] = fetchNext(f_opt); % <- This does not wait for all f_opt(i) to run! It processes them as they become available.
    
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
    f_i = D_i.fk(end);

    fs(p1_ind, p2_ind, p3_ind) = f_i;
    es(p1_ind, p2_ind, p3_ind) = D_i.e_end;
    ts(p1_ind, p2_ind, p3_ind) = D_i.t;    

    if f_i < f_best % If better performance is obtained keep this model
        D_best = D_i;
        f_best = f_i;
        ind_best = [p1_ind, p2_ind, p3_ind]; % Save indexes corresponding to best parameters
    end
end

% Prepare data for plotting - compute best value of gamma at each (t, rho) pair
gamma_best = zeros(n, n); % Store best value of gamma

for i=1:n % Parse all t values
    for k=1:n % Parse all rho values
        [~, idx] = min(fs(i, :, k));
        gamma_best(i, k) = p2(idx); % Find best gamma value
    end
end

set_plotting_parameters(1, 0); % Set text interpreter to latex

% Plot loss vs t & rho with fixed gamma -> Problem sensitivity to step-sizes
figure('Name', 'Loss vs step-size: Fixed gamma');
surf(p1, p3, squeeze(fs(:, ind_best(2), :)));
hold on;
title("\textbf{Loss vs Step-size \& Relaxation Parameter - $\gamma$ fixed}");
xlabel("$t$ - Step Size");
ylabel("$\rho$ - Relaxation Parameter");
zlabel("$f(x^{(k)}$");
grid on;
colorbar;
hold off;

% Plot best gamma vs t & rho -> gamma sensitivity to step-size
figure('Name', 'Best gamma vs step-size');
surf(p1, p3, gamma_best);
hold on;
title("\textbf{Best $\gamma$ vs Step-size \& Relaxation Parameter}");
xlabel("$t$ - Step Size");
ylabel("$\rho$ - Relaxation Parameter");
zlabel("$\gamma_{best}$");
grid on;
colorbar;
hold off;

%% Print best image
imopt_display(D_best, 'Error Evolution');
imopt_display(D_best, 'Loss Evolution');
imopt_display(D_best, 'Convergence');
imopt_display(D_best, 'Loss Convergence');

if p.save_iters % If iterates were saved
    imopt_display(D_best, 'Image Iterates', 1);
end

figure('Name', 'Deblurred Image');
imshow(D_best.xf, []);

%% Save Outputs
sweep_out = struct();
sweep_out.p1 = p1;
sweep_out.p2 = p2;
sweep_out.p3 = p3;
sweep_out.fs = fs;
sweep_out.es = es;
sweep_out.ts = ts;
sweep_out.D = D_best;

save(fn_out, "sweep_out");

%% Delete pool
delete(gcp('nocreate'));