% deblur_movie.m
%
% Makes and outputs a movie of the deblurring progression versus iteration.
%
% Inputs:
%   D: The output structure of an imopt optimization.
%
% Outputs:
%   M: The movie structure of the deblurring progress versus iteration. [struct]
%
% Author: Aidan Gerkis
% Date: 01-04-2024

function M = deblur_movie(D)
    reset(groot); % Reset ploting parameters
    % Set plotting parameters
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter', 'latex');

    % Extract vectors to plot
    xk = D.xk;
    ek = D.ek;
    e_t = D.e_end;
    
    iterations = linspace(1, length(ek), length(ek));

    % Make plots
    for i=1:length(ek) % For each iteration
        F = figure('visible', 'off');
        tiledlayout(3, 2);

        nexttile(1, [2, 2]); % Plot image - Span a 2x2 grid
        imshow(xk(:, :, i), []);
        title("\textbf{Image at Iteration " + num2str(i) + "}");
        
        nexttile(5, [1, 2]); % Plot error - Span a 1x2 grid
        plot(iterations(1:i), ek(1:i), 'LineWidth', 2);
        xlim([0, iterations(end)]);
        yline(e_t, '--r', 'LineWidth', 2)
        title("\textbf{Error Evolution vs. Iterations}");
        xlabel("Iteration");
        ylabel("Error");
        legend('', 'Error Threshold');
        grid on;

        M(i) = getframe(F); % Capture figure as a movie frame
    end

    reset(groot); % Reset ploting parameters
end    