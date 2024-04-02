% parallel_sweep
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
k = 3; % The number of parameters to sweep, each algorithm has 3 parameters

lb_1 = 0; % Lower bound on parameter one
ub_1 = 0; % Upper bound on parameter one

lb_2 = 0; % Lower bound on parameter two
ub_2 = 0; % Upper bound on parameter two

lb_3 = 0; % Lower bound on parameter three
ub_3 = 0; % Upper bound on parameter three

% Parallelization
n_pool = 14; % Number of processes to allow

% Algorithm
alg = 'primal_dr'; % Algorithm name
p_def = struct(); % Basic parameter structure to base sweep on
p_def.regularization = 'L1';
p_def.verbose = 0;
p_def.default = 0;
p_def.save_iters = 0;

%% Blur image
I = impopt_scale('cameraman.jpg'); % Import image, convert to B&W with pixels in [0,1]

kernel = fspecial('gaussian', [15,15], 5); % Define blur
b = imfilter(I,kernel); % Apply blur

noiseDensity = 0.1; % Define noise density
b = imnoise(b, 'salt & pepper', noiseDensity); % Apply noise

%% Define params to sweep
params = cell(1, n^3);

%% Run sweep
myPool = parpool(n_pool); % Initialize parallel environments

for i=1:length(params) % Initialize parallel calls
    f_opt(i) = parfeval(myPool, @imopt, 3, b, kernel, alg, params{i}); % <- These start running immediately on parfeval call
end

%% Analyze Outputs
e_best = NaN; % Set initial best error to be high

for i=1:length(params) % Loop through all iterates
    % Fetch outputs, note that the first output is an index, appended by the parfeval function
    [~, ~, e_i, D_i] = fetchNext(f_opt); % <- This does not wait for all f_opt(i) to run! It processes them as they become available.

    if e_i < e_best % If better performance is obtained keep this model
        D_best = D_i;
    end
end

%% Print best image
imopt_display(D_best, 'Error Evolution');
imopt_display(D_best, 'Convergence');
imopt_display(D_best, 'Image Iterates', 1);

figure('Name', 'Deblurred Image');
imshow(D_best.xf, []);