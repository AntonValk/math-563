function conv_rate(D)
    % conv_rate.m
    %
    % Plots convergence rate of algorithm, based on the error metric, versus iterations.
    %
    % Inputs:
    %   D: The output structure of an imopt optimization. [struct]
    %
    % Author: Aidan Gerkis
    % Date: 01-04-2024
    
    % Extract vectors to plot
    e = D.ek;
    e_t = D.e_end;

    % Set plotting parameters
    set_plotting_parameters(1, 1); % Set text and line width
    
    iterations = linspace(1, length(e), length(e));

    % Make plot
    figure('Name', 'Convergence Rates - Error');
    tiledlayout(1, 2);
    
    nexttile; % Plot semi-log of error
    semilogy(iterations, abs(e));
    hold on;
    yline(e_t, '--r')
    title("\textbf{log($\epsilon$) vs. Iterations}");
    xlabel("Iteration");
    ylabel("log($\epsilon$)");
    legend('', 'Error Threshold');
    grid on;
    hold off;

    nexttile; % Plot log-log of error
    loglog(iterations, abs(e));
    hold on;
    yline(e_t, '--r')
    title("\textbf{log($\epsilon$) vs. log(Iterations)}");
    xlabel("log(Iteration)");
    ylabel("log($\epsilon$)");
    legend('', 'Error Threshold');
    grid on;
    hold off;

    reset(groot); % Reset plotting parameters
end    