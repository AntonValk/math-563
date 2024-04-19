function stop = early_stop(res, prev, tol)
    % early_stop.m
    %
    % Adds early stop functionality to stop image algorithms iterations 
    % once progress stagnates. Implements code provided by Courtney Paquette in [1].
    %
    % Inputs:
    %   res: current residual (vector)
    %   prev: previous residual (vector)
    %
    % Outputs:
    %   stop: Binary stop / continue iterations decision (Boolean)
    %
    % Author(s): Antonios Valkanas
    % Date: 04-04-2024
    %
    % References:
    %   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
    %        in MATH 564 - Honours Convex Optimization.
    
    if norm(res - prev, 2) < tol % Check condition based on input tolerance
        
        stop = true;
    else
        stop = false;
    end
end