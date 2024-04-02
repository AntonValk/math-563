% benchmark_parallel.m
%
% Analyzes a functions parallel computing performance as a function of the
% number of pools. Considers compute time, speed up, and data transfer.
%
% Author: Aidan Gerkis
% Date: 29-02-2024

clear; clc; close all;

%% Parameters
% The function to analyze, must take the number of pools as an input and 
% output sim_time [scalar], speed_up ratio [scalar], and bytes transferred
% [1 x 2 array]
fn_handle = @(n_pool)parallel_sweep_benchmark(n_pool);

% Pool Sizes to Consider
n_pool = [2, 4, 8, 12, 16, 20, 24];

% Initialize Output Arrays
sim_times = zeros(1, length(n_pool));
speed_ups = zeros(1, length(n_pool));
bytes_tf = zeros(2, length(n_pool));

%% Simulate With Different n_pool
for i=1:length(n_pool)
    [st, su, b_tf] = fn_handle(n_pool(i)); % Simulate
    
    disp("Simulation of function with " + num2str(n_pool(i)) + " workers completed in " + num2str(st) + "s.");
    disp("");

    % Save Performance Metrics
    sim_times(i) = st;
    speed_ups(i) = su;
    bytes_tf(:, i) = b_tf';
end

eff = (speed_ups./n_pool)*100; % Compute efficiency as the proportion of speed up to number of workers (ideal efficiency: each worker takes exactly the same amount of time to compute)

%% Plots
close all; % Close any plots made by the function being benchmarked
set_plotting_parameters(1, 1);

figure; % Plot speed up ratio & efficiency vs. pool size
plot(n_pool, speed_ups, '-*');
hold on;
title("\textbf{Speed-Up \& Efficiency vs. Number of Workers}");
xlabel("Number of Workers");
ylabel("Speed-Up Ratio");
yyaxis right;
plot(n_pool, eff, '--o');
ylabel("Efficiency [\%]");
grid on;
hold off;

figure; % Plot sim time vs. pool size
plot(n_pool, sim_times);
hold on;
title("\textbf{Simulation Time vs. Number of Workers}");
xlabel("Number of Workers");
ylabel("Simulation Time [s]");
grid on;
hold off;

figure; % Plot transfer size vs. pool size
plot(n_pool, bytes_tf(1, :), '-*');
hold on;
title("\textbf{Data Transferred vs. Number of Workers}");
xlabel("Number of Workers");
ylabel("Client$\rightarrow$Workers [GB]");
yyaxis right;
plot(n_pool, bytes_tf(2, :), '--o');
ylabel("Workers$\rightarrow$Client [GB]");
grid on;
hold off;