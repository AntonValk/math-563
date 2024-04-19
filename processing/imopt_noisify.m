function n = imopt_noisify(im, type, varargin)
    % imopt_noisify.m
    %
    % Apply noise to an image. Uses salt & pepper with density 0.1 as a
    % default.
    %
    % Inputs:
    %   im: The image to apply noise too. [m x n Matrix]
    %   type: The type of noise to apply. May be one of:
    %               - 'gaussian'
    %               - 'localvar'
    %               - 'poisson'
    %               - 'salt & pepper'
    %               - 'speckle'
    %   varargin: Contains the parameters for the specified noise type, see
    %             'imnoise' documentation for details. Contains at least 0
    %             values and at most 2.
    %
    % Outputs:
    %   n: The image with noise applied. [m x n Matrix]
    %
    % Author(s): Aidan Gerkis
    % Date: 04-04-2024
    
    if ~isa(im, 'numeric') % Validate image input
        disp("Error in function call, expected a numeric but got a " + class(im));
    end
    
    if nargin >= 2 % Validate type input, if it was passed
        if ~isa(type, 'char')
            disp("Error in function call, expected a char array but got a " + class(type));
        end
    end
    
    switch length(varargin) % Build function handle for denoising
        case 0 % Either 1 or 2 inputs were 
            if nargin == 1 % If no noise type was specified then use the default type
                noisify = @(im) imnoise(im, 'salt & pepper', 0.1);
            else
                noisify = @(im) imnoise(im, type);
            end
        case 1
            noisify = @(im) imnoise(im, type, varargin{1});
        case 2
            noisify = @(im) imnoise(im, type, varargin{1}, varargin{2});
    end

    n = noisify(im); % Apply noise
end



