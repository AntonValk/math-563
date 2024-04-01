% imopt.m
%
% THE package.
%
% Inputs:
%   Some things.
%
% Outputs:
%   some (hopefully less blurry) things
%
% Author(s): Aidan Gerkis
% Date: 30-3-2024

function varargout = imopt(b, kernel, alg, p_in)
    in_names = ["b", "kernel", "alg", "p_in"]; % List of names of inputs
    input_types = ["numeric", "numeric", "char", "struct"]; % List of allowed types for inputs

    fields_allowed = ["verbose", "display", "save_iters", "max_iter", "e_t", ...
        "e_meas", "regularization", "t", "s", "gamma", "rho"]; % List of all legal field names for the parameter structure
    fields_type = ["logical", "logical", "logical", "numeric", "numeric", ...
        "char", "char", "numeric", "numeric", "numeric", "numeric"]; % List of types for all fields
    
    %% TODO: Add input parsing -> Make sure inputs are of the correct type
    
    % Parse Inputs - Check number & assign default values
    switch nargin
        case 1
            disp("Error in imopt: Error in function call, too few input arguments.");
        case 2
            argin = cell(1, nargin);
            argin{1} = b;
            argin{2} = kernel;
            alg = 'primal_dr'; % Assign default algorithm
        case 3
            argin = cell(1, nargin);
            argin{1} = b;
            argin{2} = kernel;
            argin{3} = alg;
        case 4
            argin = cell(1, nargin);
            argin{1} = b;
            argin{2} = kernel;
            argin{3} = alg;
            argin{4} = p_in;
        otherwise
            disp("Error in imopt: Error in function call, too many input arguments.");
    end

    % Parse inputs - check type
    for i=1:nargin
        if ~isa(argin{i}, input_types(i)) % Display error message and exit if type is incorrect
            disp("Error in imopt: Unexpected type for parameter '" + in_names(i) + ...
                "'. Expected a " + input_types(i) + " but got a " +class(argin{i}) + ".");
            return
        end
    end

    % Evaluate input parameter structure
    params = get_default(alg, b);

    if exist('p_in') % If custom parameters were passed overwrite default values
        names = fieldnames(p_in);

        for i=1:length(names) % Overwrite defaults if a legal field was passed
            if ismember(names{i}, fields_allowed) % Check if parameter name is correct
                if isa(p_in.(names{i}), fields_type(fields_allowed == names{i})) % Check if parameter type is correct
                    params.(names{i}) = p_in.(names{i});
                else
                    disp("Error in imopt: Unexpected type for parameter '" + names{i} + ...
                        "'. Expected a " + fields_type(fields_allowed == names{i}) + ...
                        " but got a " + class(p_in.(names{i})));
                    return
                end
            else
                disp("Error in imopt: Unrecognized parameter passed in p_in.")
                return
            end
        end
    end
  
    
    % NOTE: SWAP ROLES OF s and t in Chambolle Pock! For consistency with
    % other algs....
    
    %% TODO: How to handle different error metrics....

    if params.verbose
        disp("==============Parsing Input Parameters==============");
    end

    % Build function handle for prox g
    prox_g = @(x, t) prox_tg(x, t, params.gamma, params.regularization, b);
    
    %% TODO: Double check which variables are actually needed in the function
    %%       Address mycourses posts comment -> Should be easy with this implementation

    % Build function handle for algorithm & structure to store inputs    
    switch alg
        case 'primal_dr'
            deblur = @(im)primal_douglasrachford_splitting(im, kernel, params.x_init, prox_g, params.t, params.rho, params.max_iter, params.e_t, params.save_iters);
            
            % Compile parameters for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "rho: " + num2str(params.rho);
            alg_name = "Primal Douglas-Rachford Splitting";
        case 'primaldual_dr'
            deblur = @(im)primaldual_douglasrachford_splitting(im, kernel, params.x_init, params.t, params.rho, params.max_iter, params.e_t, params.save_iters);
            
            % Compile parameters for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "rho: " + num2str(params.rho);
            alg_name = "Primal-Dual Douglas-Rachford Splitting";
        case 'admm'
            %% TODO: Make sure params are right in this function call
            deblur = @(im)admm(im, kernel, params.x_init, params.rho, params.max_iter, params.e_t, params.save_iters);
            
            % Compile outputs for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "rho: " + num2str(params.rho);
            alg_name = "Alternating Direction Method of Multipliers";
        case 'chambolle_pock'
            deblur = @(im)chambolle_pock(im, kernel, params.x_init, params.t, params.s, params.max_iter, params.e_t, params.save_iters);
            
            % Compile outputs for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "s: " + num2str(params.s);
            alg_name = "Chambolle Pock";
        otherwise %% Change to pass error
            disp("Error in imopt: Unknown algorithm specified.");
            return
    end

    % Run deblurring algorithm
    disp("==============Running Deblurring Algorithm==============");

    if params.verbose
        disp("Algorithm: " + alg_name);
        disp("  " + p_vals(1));
        disp("  " + p_vals(2));
        disp("  Regularization: " + params.regularization);
        disp("  Max Iterations: " + params.max_iter);
        disp("  Error Threshold: " + params.e_t);
    end

    D = deblur(b);

    disp("==============Deblurring Completed==============");

    % Print outputs
    if params.verbose
        disp("Deblurring algorithm completed successfully.");
        disp("Total Time: " + num2str(D.t) + " s");
        disp("Total Iterations: " + num2str(D.k_end));
        disp("Final Error: " + num2str(D.e_end));

        if rand >= 0.05
            disp("IMOPT - Math 563 Final Project");
        else
            disp("Thanks for using IMOPT, have a great day!");
        end
    end
    
    % Display output plots
    if params.display
        imopt_display(D, 'Error Evolution');
        imopt_display(D, 'Convergence');
        imopt_display(D, 'Image Iterates', 1);
    end
    
    % Compile outputs
    D.inputs = params;
    D.inputs.alg_name = alg_name;

    varargout = cell(1, nargout);

    switch nargout % Create outputs and return
        case 1 % Return only the final image
            varargout{1} = D.xf;
        case 2 % Return final image and final error
            varargout{1} = D.xf;
            varargout{2} = D.ef;
        case 3 % Return final image, final error, and output structure
            varargout{1} = D.xf;
            varargout{2} = D.e_end;
            varargout{3} = D;
        otherwise
            disp("Error in imopt: Too many outputs requested.");
    end
end