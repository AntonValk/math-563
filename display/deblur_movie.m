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
    % Set plotting parameters
    set_plotting_parameters(1, 0); % Only set text to latex

    % Extract vectors to plot
    xk = D.xk;
    fk = D.fk;
    
    iterations = linspace(1, length(fk), length(fk));

    % Make plots
    for i=1:length(fk) % For each iteration
        F = figure('visible', 'off');
        tiledlayout(3, 2);

        nexttile(1, [2, 2]); % Plot image - Span a 2x2 grid
        imshow(xk(:, :, i), []);
        title("\textbf{Image at Iteration " + num2str(i) + "}");
        
        nexttile(5, [1, 2]); % Plot error - Span a 1x2 grid
        semilogy(iterations(1:i), fk(1:i), 'LineWidth', 2);
        xlim([0, iterations(end)]);
        ylim([1, max(fk)]);
        title("\textbf{Loss vs. Iteration}");
        xlabel("Iteration");
        ylabel("Loss");
        grid on;

        M(i) = getframe(F); % Capture figure as a movie frame
    end

    reset(groot); % Reset ploting parameters
end    