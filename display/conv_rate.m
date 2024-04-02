% conv_rate.m
%
% Plots convergence rate of algorithm versus iterations.
%
% Inputs:
%   D: The output structure of an imopt optimization. [struct]
%
% Author: Aidan Gerkis
% Date: 01-04-2024

function conv_rate(D)
    % Extract vectors to plot
    e = D.ek;
    e_t = D.e_end;

    % Set plotting parameters
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot,'DefaultLineLineWidth',2);
    set(groot,'DefaultContourLineWidth',1.5);
    set(groot,'DefaultFunctionContourLineWidth',2);
    
    iterations = linspace(1, length(e), length(e));

    % Make plot
    figure('Name', 'Convergence Rates');
    tiledlayout(1, 2);
    
    nexttile; % Plot semi-log of error
    semilogy(iterations, e);
    hold on;
    yline(e_t, '--r')
    title("\textbf{log($\epsilon$) vs. Iterations}");
    xlabel("Iteration");
    ylabel("log($\epsilon$)");
    legend('', 'Error Threshold');
    grid on;
    hold off;

    nexttile; % Plot log-log of error
    loglog(iterations, e);
    hold on;
    yline(e_t, '--r')
    title("\textbf{log(Error Evolution) vs. log(Iterations)}");
    xlabel("log(Iteration)");
    ylabel("log(Error)");
    legend('', 'Error Threshold');
    grid on;
    hold off;

    reset(groot); % Reset ploting parameters
end    