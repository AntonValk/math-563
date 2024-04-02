% error_evo.m
%
% Plots error evolution versus iterations.
%
% Inputs:
%   D: The output structure of an imopt optimization. [struct]
%
% Author: Aidan Gerkis
% Date: 01-04-2024

function error_evo(D)
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

    % Make plot
    iterations = linspace(1, length(e), length(e));
    figure('Name', 'Error Evolution');
    plot(iterations, e)
    hold on;
    yline(e_t, '--r')
    title("\textbf{$\epsilon$ vs. Iterations}");
    xlabel("Iteration");
    ylabel("$\epsilon$");
    legend('', 'Error Threshold');
    grid on;
    hold off;

    reset(groot); % Reset ploting parameters
end    