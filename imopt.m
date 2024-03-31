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

function [xf, ef, D] = imopt(b, kernel, alg, p_in)
    switch nargin % Parse input algorithm
        case 2 % Use default algorithm
            alg = "default alg";
        otherwise % Do nothing
    end

    % Evaluate input parameter structure
    params = get_default();

    if exists(p_in) % If custom parameters were passed overwrite default values
        names = fieldnames(p_in);

        for i=1:length(names)
            params.(names{i}) = p_in.(names{i});
        end
    end
    
    % Add input parsing -> Make sure inputs are not incorrect
    
    % NOTE: SWAP ROLES OF s and t in Chambolle Pock! For consistency with
    % other algs....
    
    %% TODO: How to handle different error metrics....
    
    v = params.verbose;
    d = params.display;

    if v
        disp("==============Parsing Input Parameters==============");
    end

    % Build function handle for prox g
    prox_g = @(x) prox_tg(x, params.t, params.gamma, params.regularization, b);
    
    %% TODO: Double check which variables are actually needed in the function
    %%       Address mycourses posts comment -> Should be easy with this implementation
    %%       Change function handles of the algs

    % Build function handle for algorithm
    switch alg
        case 'primal_dr'
            deblur = @(im)primal_douglasrachford_splitting(im, kernel, params.rho, params.max_iter, params.err_thresh);
            
            % Compile outputs for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "rho: " + num2str(params.rho);
            alg_name = "Primal Douglas-Rachford Splitting";
        case 'primaldual_dr'
            deblur = @(im)primaldual_douglasrachford_splitting(im, kernel, params.rho, params.max_iter, params.err_thresh);
            
            % Compile outputs for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "rho: " + num2str(params.rho);
            alg_name = "Primal-Dual Douglas-Rachford Splitting";
        case 'admm'
            %% TODO: Make sure params are right in this function call
            deblur = @(im)admm(im, kernel, params.rho, params.max_iter, params.err_thresh);
            
            % Compile outputs for verbose mode
            p_vals(1) = "t: " + num2str(params.t);
            p_vals(2) = "rho: " + num2str(params.rho);
            alg_name = "Alternating Direction Method of Multipliers";
        case 'chambolle_pock'
            deblur = @(im)chambolle_pock(im, kernel, params.t, params.s, params.max_iter, params.err_thresh);
            
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

    if v
        disp("Algorithm: " + alg_name);
        disp("  " + p_vals(1));
        disp("  " + p_vals(2));
        disp("  Regularization: " + params.regularization);
        disp("  Max Iterations: " + params.max_iter);
        disp("  Error Threshold: " + params.err_thresh);
    end

    D = deblur(b);
    disp("==============Deblurring Completed==============");

    % Print outputs
    if v
        if D.success
            disp("Deblurring algorithm completed successfully.");
        else
            disp("Deblurring algorithm failed.");
        end

        disp("Total Time: " + num2str(D.t) + " s");
        disp("Total Iterations: " + num2str(D.k_end));
        disp("Final Error: " + num2str(D.e_end));

        if rand >= 0.05
            disp("IMOPT - Math 563 Final Project");
        else
            disp("Thanks for using IMOPT, have a great day!");
        end
    end
    
    % Display output statistics
    if d
        %% TODO: Display code
        % Write these as functions that plot the desired value
        % Error evolution
        % Convergence rate (linear & sublinear?)
        % Movie of image evolution over time (if save is on)
    end

    switch nargout % Create outputs and return
        case 1 % Return only the final image
            xf = D.xf;
        case 2 % Return final image and final error
            xf = D.xf;
            ef = D.ef;
        case 3 % Return final image, final error, and output structure
            xf = D.xf;
            ef = D.ef;
        otherwise
            disp("Error in imopt: Too many outputs requested.");
    end
end