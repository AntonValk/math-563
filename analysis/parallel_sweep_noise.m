% parallel_sweep_noise.m
%
% Performs a sweep of the three model hyperparameters, utilizing MATLAB's
% parallelization toolbox. Requires that the parallel computing toolbox be
% installed.
%
% Author: Aidan Gerkis
% Date:02-04-2024

clear; clc; close all;

%% Parameters
n = 100; % Number of parameters to sweep in each direction
k = 2; % The number of parameters to sweep, each algorithm has 3 parameters

lb_1 = 1e-15; % Lower bound on parameter one (noise density)
ub_1 = 0.85; % Upper bound on parameter one (noise density)

lb_2 = 1e-15; % Lower bound on parameter two (gamma)
ub_2 = 5; % Upper bound on parameter two (gamma)

fn_out = "n10-noisesweep-primal";

% Parallelization
%n_pool = 16; % Number of processes to allow. Chosen from benchmarking on a 24 core computer
n_pool = 8; % Number of processes to allow. Use this for your PCs! Or maybe even lower, check your system info.

% Algorithm
alg = 'primal_dr'; % Algorithm name
p = struct(); % Basic parameter structure to base sweep on
p.regularization = 'L1';
p.verbose = false;
p.display = false;
p.save_iters = false;
p.tol = 0.1; % Enable early stop, for the love of god, enable early stop

%% Blur image
kernel = fspecial('gaussian', [9, 9], 4); % Define kernel
I_true = imopt_scale('cameraman.jpg'); % Import image, convert to B&W with pixels in [0,1]
blurred = imopt_blur(I_true, kernel);

%% Define params to sweep
p1 = linspace(lb_1, ub_1, n); % Values for parameter 1
p2 = linspace(lb_2, ub_2, n); % Values for parameter 2

%% Run sweep
myPool = parpool(n_pool); % Initialize parallel environments
idx = 1; % Index parallel calls

% Initialize parallel calls, loop backwards so all memory is allocated at
% first step
for i=n:-1:1 % Sweep values of parameter 1
    for j=n:-1:1 % Sweep values of parameter 2
            % Assign parameter values
            p.gamma = p2(j);

            % Blur image
            b = imopt_noisify(blurred, 'salt & pepper', p1(i));
            b = cat(3, b, I_true);

            f_opt(idx) = parfeval(myPool, @imopt, 3, b, kernel, alg, p); % <- These start running immediately on parfeval call
            idx = idx + 1;
    end
end

%% Analyze Outputs
fs = zeros(n, n); % Store loss for each parameter combination
es = zeros(n, n); % Store error for each parameter combination
ts = zeros(n, n); % Store computation time for each parameter combination

f_best = Inf; % Set initial best loss to be high

for i=1:length(f_opt) % Loop through all iterates
    % Fetch outputs, note that the first output is an index, appended by the parfeval function
    [idx, x, e_i, D_i] = fetchNext(f_opt); % <- This does not wait for all f_opt(i) to run! It processes them as they become available.
    
    % Extract parameters of current run
    p1_ind = floor((idx - 1)/10) + 1;
    p2_ind = find(p2 == D_i.inputs.gamma);
    
    % Save performance metrics
    f_i = D_i.fk(end);

    fs(p1_ind, p2_ind) = f_i;
    es(p1_ind, p2_ind) = D_i.e_end;
    ts(p1_ind, p2_ind) = D_i.t;    

    if f_i < f_best % If better performance is obtained keep this model
        D_best = D_i;
        f_best = f_i;
        ind_best = [p1_ind, p2_ind]; % Save indexes corresponding to best parameters
    end
end

% Prepare data for plotting - compute best value of gamma at each noise
% density value
gamma_best = zeros(n); % Store best value of gamma

for i=1:n % Parse all noise density values
    [~, idx] = min(fs(i, :));
    gamma_best(i) = p2(idx); % Find best gamma value
end

set_plotting_parameters(1, 0); % Set text interpreter to latex

% Plot loss vs gamma and noise density -> Problem sensitivity to gamma as
% noise density changes
figure('Name', 'Loss heatmap'); % As heatmap
heatmap(p1, p2, squeeze(fs), 'colormap', parula, 'celllabelcolor', 'none');
title("Loss vs noise-density & \gamma");
xlabel("t");
ylabel("\rho");

figure('Name', 'Error heatmap'); % Error as heatmap
heatmap(p1, p2, squeeze(es), 'colormap', parula, 'celllabelcolor', 'none');
title("Error vs t & \rho");
xlabel("t");
ylabel("\rho");

% Plot best gamma vs t & rho -> gamma sensitivity to step-size
figure('Name', 'Gamma surface');
plot(p1, gamma_best);
hold on;
title("\textbf{Best $\gamma$ vs Noise Density}");
xlabel("Noise Density");
ylabel("$\gamma$");
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