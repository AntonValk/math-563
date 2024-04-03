% loss_conv_rate.m
%
% Plots convergence rate of algorithm, based on the loss, versus iterations.
% Requires that the optimal point is given in D.inputs. If it is not then
% D.xf is used.
%
% Inputs:
%   D: The output structure of an imopt optimization. [struct]
%
% Author: Aidan Gerkis
% Date: 03-04-2024

function loss_conv_rate(D)
    % Extract vectors to plot
    f = D.fk;

    if isfield(D.inputs, 'x_true') % If true image was passed
        fs = imopt_loss(D.inputs.x_true, D.inputs.b, D.inputs.gamma, D.inputs.regularization, D.inputs.kernel);
    else % Use final output as the "true" image
        fs = imopt_loss(D.xf, D.inputs.b, D.inputs.gamma, D.inputs.regularization, D.inputs.kernel);
    end

    % Set plotting parameters
    set_plotting_parameters(1, 1); % Set text and line width
    
    iterations = linspace(1, length(f), length(f));

    % Make plot
    figure('Name', 'Convergence Rates - Loss');
    tiledlayout(1, 2);
    
    nexttile; % Plot semi-log of loss - fs
    semilogy(iterations, abs(f - fs));
    hold on;
    title("\textbf{log($f(x^{(k)}) - f^*$) vs. Iterations}");
    xlabel("Iteration");
    ylabel("log($f(x^{(k)}) - f^*$)");
    grid on;
    hold off;

    nexttile; % Plot log-log of error
    loglog(iterations, f - fs);
    hold on;
    title("\textbf{log($f(x^{(k)}) - f^*$) vs. log(Iterations)}");
    xlabel("log(Iteration)");
    ylabel("log($f(x^{(k)}) - f^*$)");
    grid on;
    hold off;

    reset(groot); % Reset plotting parameters
end    