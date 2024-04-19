function loss_evo(D)
    % loss_evo.m
    %
    % Plots loss versus iterations.
    %
    % Inputs:
    %   D: The output structure of an imopt optimization. [struct]
    %
    % Author: Aidan Gerkis
    % Date: 03-04-2024
    % Extract vectors to plot
    f = D.fk;

    % Set plotting parameters
    set_plotting_parameters(1, 1); % Set text and line width

    % Make plot
    iterations = linspace(1, length(f), length(f));
    figure('Name', 'Loss Evolution');
    plot(iterations, f)
    hold on;
    title("\textbf{$f(x^{(k)}$ vs. Iterations}");
    xlabel("Iteration");
    ylabel("$f(x^{(k)}$");
    grid on;
    hold off;

    reset(groot); % Reset plotting parameters
end    