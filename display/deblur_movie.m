function M = deblur_movie(D)
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
    
    % Set plotting parameters
    set_plotting_parameters(1, 0); % Only set text to latex

    % Extract vectors to plot
    xk = D.xk;
    fk = D.fk;
    ns = 1;

    if D.inputs.save_iters == 2 % Trim loss vector if sparse save was used
        N = linspace(1, length(fk), length(fk));
        fk = fk(mod(N-1, D.inputs.ns) == 0);
        ns = D.inputs.ns;
    end
    
    iterations = linspace(1, length(fk), length(fk));
    
    figure('visible', 'off'); % Get y-limits to use when plotting loss vs. iterations
    loglog(iterations(1:length(fk)).*ns, fk, 'LineWidth', 2);
    yl = ylim;

    % Make plots
    for i=length(fk):-1:1 % For each iteration. Loop backwards for better memory allocation
        F = figure('visible', 'off');
        tiledlayout(3, 2);

        nexttile(1, [2, 2]); % Plot image - Span a 2x2 grid
        imshow(xk(:, :, i), []);
        title("\textbf{Image at Iteration " + num2str(i*ns) + "}");
        
        nexttile(5, [1, 2]); % Plot error - Span a 1x2 grid
        loglog(iterations(1:i).*ns, fk(1:i), 'LineWidth', 2);
        xlim([0, iterations(length(fk)).*ns]);
        ylim(yl);
        title("\textbf{Loss vs. Iteration}");
        xlabel("Iteration");
        ylabel("Loss");
        grid on;

        M(i) = getframe(F); % Capture figure as a movie frame
    end
    
    flip(M); % Need to flip M as it was constructed in reverse order
    reset(groot); % Reset ploting parameters
end    