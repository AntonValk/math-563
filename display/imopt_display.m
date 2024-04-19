function imopt_display(D, plot, n)
    % imopt_display.m
    %
    % A wrapper function that makes the requested plot type, given the output
    % structure of an imopt optimization.
    %
    % Inputs:
    %   D: The output structure of an imopt optimization. [struct]
    %   plot: The type of plot to display. [char array]
    %   n: The number of times to play the movie. [integer]
    %   
    % Author: Aidan Gerkis
    % Date: 01-04-2024
    
    switch plot
        case 'Error Evolution' % Make plot of error evolution
            error_evo(D);
        case 'Convergence' % Make plot of convergence
            conv_rate(D);
        case 'Loss Evolution' % Plot loss vs iterates
            loss_evo(D);
        case 'Loss Convergence' % Make plots of convergence wrt loss
            loss_conv_rate(D);
        case 'Image Iterates' % Make and play movie
            M = deblur_movie(D); % Make Movie

            % Play movie
            [h, w, ~] = size(M(1).cdata);
            f = figure('Name', 'Image Evolution vs. Iteration');
            set(f, 'position', [0 0 1.25*w 1.25*h]);
            axis off;
            movie(M, n);
        otherwise
            error("Unknown plot type requested");
    end
end