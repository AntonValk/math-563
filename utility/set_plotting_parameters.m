function set_plotting_parameters(txt, plt)
    % set_plotting_parameters.m
    %
    % Sets my desired graphical display parameters for generating MATLAB plots.
    % First resets the parameters to default.
    %
    % Inputs:
    %   txt: A boolean indicating whether to set text parameters (0 - Don't set, 1 - Set)
    %   plt: A boolean indicating whether to set plot parameters (0 - Don't set, 1 - Set)
    
    reset(groot);
    
    if txt % Change Text Interpreter
        set(groot, 'defaultAxesTickLabelInterpreter','latex');
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaultTextInterpreter', 'latex');
    end
    
    if plt % Change Line Style
        set(groot,'DefaultLineLineWidth',2);
        set(groot,'DefaultContourLineWidth',1.5);
        set(groot,'DefaultFunctionContourLineWidth',2);
    end
end