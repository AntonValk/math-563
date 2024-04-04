% imopt.m
%
% A package to implement four image deblurring algorithms: primal
% douglas-rachford splitting, primal-dual douglas-rachford splitting, admm,
% and chambolle-pock. Supports a very of interactability, allowing the user
% access to all algorithm hyperparameters. Also supports different outputs
% and display modes.
%
% May be called in one of three ways:
%   imopt(b, kernel)
%   imopt(b, kernel, alg)
%   imopt(b, kernel, alg, p_in)
%
% And supports up to three outputs.
%
% Inputs:
%   b: The blurred & true images. Blurred in first index, true in second. [m x n x 2 Tensor]
%   kernel: The kernel used to deblur the image. [k x k Matrix]
%   alg: The algorithm with which to deblur the image. Should be one of:
%           - 'primal_dr'
%           - 'primaldual_dr'
%           - 'admm'
%           - 'chambolle_pock'
%   p_in: A parameter structure containing inputs to the algorithms. May
%   contain any (or none) of the following:
%           x_init: The initial guess for the deblurred image. [m x n Matrix]
%           t: Step size. [Double] (For all algorithms)
%           s: Step size. [Double] (For the Chambolle-Pock method)
%           g: The constant modifying the iso-norm in the problem statement. [Double]
%              (For all algorithms)
%           rho: Regularization parameter. [Double] (For the primal-dr,
%                primaldual-dr, and admm methods.)
%           regularization: The regularization type to use. May be one of
%                            - 'L1'
%                            - 'L2'
%           metric: The error metric to use. May be one of
%                            - 'rmse'
%                            - 'psnr'
%                            - 'varinfo'
%           L: The number of bins for the variation of information metric. [Double]
%              Mandatory if metric = 'varinfo'.
%           max_iters: Maximum number of iterations. [Integer]
%           e_t: Error threshold. [Double]
%           verbose: A boolean, indicating whether verbose outputs should be printed. [Logical]
%           display: A boolean, indicating whether plots should be made. [Logical]
%           save_iters: A boolean, indicating whether the image iterates should be saved. [Logical]
%
% Outputs:
%   xf: The final image, after deblurring. [m x n Matrix]
%   e_end: The loss function value at the final iteration. [Double]
%   D: The complete output structure from the algorithm calls. [Struct]
%       xf: The final image. [m x n Matrix]
%       t: The time it took to run the optimization algorithm. [Double, Seconds]
%       k_end: The number of iterations ran. [Integer]
%       e_end: Error at the final iteration. [Double]
%       ek: Error at each iteration [1 x k_end Matrix]
%       xk: The image at each iteration [m x n x k_end Matrix] (Only output if save is true)
%       inputs: A structure containing the inputs to the imopt package. [Structure]
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
    
    % Parse Outputs - Check number of outputs requested
    if nargout > 3
        error("Too many outputs requested.");
    end

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
            error("Error in function call, too many input arguments.");
    end

    % Parse inputs - check type
    for i=1:nargin
        if ~isa(argin{i}, input_types(i)) % Display error message and exit if type is incorrect
            error("Unexpected type for parameter '" + in_names(i) + ...
                "'. Expected a " + input_types(i) + " but got a " +class(argin{i}) + ".");
        end
    end
    
    % Extract blurred image and ground truth
    [~, ~, k] = size(b); 
    if k == 2
        x_true = b(:, :, 2);
        b = b(:, :, 1);
    else
        error("Error in function call, expected ground truth to be passed.");
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
                    error("Unexpected type for parameter '" + names{i} + ...
                        "'. Expected a " + fields_type(fields_allowed == names{i}) + ...
                        " but got a " + class(p_in.(names{i})));
                end
            else
                error("Unrecognized parameter passed in p_in.");
            end
        end
    end
    
    if params.verbose
        disp("================Parsing Input Parameters================");
    end
    
    f_handles = struct(); % Compile all functions into structure

    % Build error metric function
    switch params.metric
        case 'rmse'
            f_handles.err_eval = @(x_hat) imopt_rmse(x_true, x_hat);
        case 'psnr'
            f_handles.err_eval = @(x_hat) imopt_psnr(x_true, x_hat);
        case 'varinfo'
            if isfield(params, 'L')
                f_handles.err_eval = @(x_hat) imopt_var_info(x_true, x_hat, params.L);
            else
                error("Error when calling 'imopt_var_info'. Number of bins, 'L', was not specified.");
            end
        otherwise
            error("Unknown metric type specified");
    end

    % Build loss function
    f_handles.f_loss = @(x)imopt_loss(x, b, params.gamma, params.regularization, kernel);

    % Build generic function handle for prox of regularization function
    f_handles.prox_l = @(x, t) prox_ln(x, t, params.regularization, b);

    % Build function handle for algorithm & structure to store inputs    
    switch alg
        case 'primal_dr'
            deblur = @(im)primal_douglasrachford_splitting(im, kernel, params.x_init, f_handles, params.t, params.gamma, params.rho, params.max_iter, params.e_t, params.save_iters, params.verbose);
            
            % Compile parameters for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "rho: " + num2str(params.rho);
            alg_name = "Primal Douglas-Rachford Splitting";
        case 'primaldual_dr'
            deblur = @(im)primaldual_douglasrachford_splitting(im, kernel, params.x_init, f_handles, params.t, params.gamma, params.rho, params.max_iter, params.e_t, params.save_iters, params.verbose);

            % Compile parameters for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "rho: " + num2str(params.rho);
            alg_name = "Primal-Dual Douglas-Rachford Splitting";
        case 'admm'
            deblur = @(im)admm(im, kernel, params.x_init, f_handles, params.t, params.rho, params.gamma, params.max_iter, params.e_t, params.save_iters, params.verbose);
            
            % Compile outputs for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "rho: " + num2str(params.rho);
            alg_name = "Alternating Direction Method of Multipliers";
        case 'chambolle_pock'
            deblur = @(im)chambolle_pock(im, kernel, params.x_init, f_handles, params.t, params.s, params.gamma, params.max_iter, params.e_t, params.save_iters, params.verbose);
            
            % Compile outputs for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "s: " + num2str(params.s);
            alg_name = "Chambolle Pock";
        otherwise %% Change to pass error
            error("Unknown optimization algorithm specified.");
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

    disp("==================Deblurring Completed==================");

    % Compile outputs
    D.inputs = params;
    D.inputs.b = b;
    if exist('x_true')
        D.inputs.x_true = x_true;
    end
    D.inputs.alg_name = alg_name;
    D.inputs.kernel = kernel;

    % Print outputs
    if params.verbose
        disp("Deblurring algorithm completed successfully.");
        disp("Total Time: " + num2str(D.t) + " s");
        disp("Total Iterations: " + num2str(D.k_end));
        disp("Final Error: " + num2str(D.e_end));
        disp("Final Loss: " + num2str(D.fk(end)));

        if rand >= 0.05
            disp("IMOPT - Math 563 Final Project");
        else
            disp("Thanks for using IMOPT, have a great day!");
        end
    end

    % Display output plots
    if params.display
        imopt_display(D, 'Error Evolution');
        imopt_display(D, 'Loss Evolution');
        imopt_display(D, 'Convergence');
        imopt_display(D, 'Loss Convergence');
        if params.save_iters
            imopt_display(D, 'Image Iterates', 1);
        end
    end
    
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
        % otherwise case is handled on function call
    end
end