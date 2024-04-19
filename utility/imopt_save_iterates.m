function imopt_save_iterates(D, dir, name)
    % imopt_save_iterates.m
    % 
    % Saves the image iterates from each algorithm iteration as pngs in the
    % specified directory and with the specified filename. Only works if the
    % image iterates were aved in the algorithm call.
    %
    % Inputs:
    %   D: The output structure of an imopt call. [Structure]
    %   dir: The directory to save the image. [String] (Optional)
    %   name: The name to save the image under.[String] (Optional)
    %
    % Author: Aidan Gerkis
    % Date: 15-04-2025
    
    input_names = ["D", "dir", "name"];
    input_types = ["struct", "string", "string"];

    switch nargin
        case 1 % Use default directory and name, or the ones specified by the user when calling imopt
            dir = D.inputs.dir;
            name = D.inputs.im_name;

            % Save arguments in cell array for parsing
            argin{1} = D;
        case 2 % Use default name, directory provided
            name = D.inputs.im_name;

            % Save arguments in cell array for parsing
            argin{1} = D;
            argin{2} = dir;
        case 3 % All parameters specified
            % Save arguments in cell array for parsing
            argin{1} = D;
            argin{2} = dir;
            argin{3} = name;
        otherwise
            error("Error in function call, too many input arguments.");
    end

    % Parse inputs - check type
    for i=1:nargin
        if ~isa(argin{i}, input_types(i)) % Display error message and exit if type is incorrect
            error("Unexpected type for parameter '" + input_names(i) + ...
                "'. Expected a " + input_types(i) + " but got a " + class(argin{i}) + ".");
        end
    end
    
    % Ensure image iterates were saved
    if D.inputs.save_iters == 0
        error("Image iterates not saved, exitting imopt_save_iterates.");
    end

    % Save iterates as pngs
    n = length(D.xk(1, 1, :)); % Get number of iterates to save
    ns = 1; % Scaling factor for iterate number

    if D.inputs.save_iters == 2 % If sparse saving was used adjust ns
        ns = D.inputs.ns;
    end

    for i=1:n % Save iterates
        iname = dir + "\" + name + "_iterate" + num2str(i*ns) + ".png";
        imwrite(squeeze(D.xk(:, :, i)), iname);
    end
end