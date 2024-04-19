function b = imopt_blur(im, kernel, padding)
    % imopt_blur.m
    %
    % Blurs an image using padding to remove rough edges that arise when
    % post-processing the image. Enforces 'convolution' as the filtering
    % method.
    %
    % Inputs:
    %   im: The image to blur. [m x n Matrix]
    %   kernel: The kernel to use when applying blur. [k x k Matrix] (Optional)
    %   padding: The padding type to use. May be one of: (Optional)
    %               - 'symmetric'
    %               - 'replicate'
    %               - 'circular'
    %
    % Outputs:
    %   b: The blurred image. [m x n Matrix]
    %
    % Usage:
    %   b = imopt_blur(im); % Blurs the input image with a 9x9 gaussian kernel
    %                         of sigma = 4. Uses replicate padding.
    %   b = imopt_blur(im, my_kernel); % Blurs the input image with kernel
    %                                    'my_kernel' and replicate padding.
    %   b = imopt_blur(im, my_kernel, my_padding); % Blurs the input image with
    %                                                kernel 'my_kernel' and
    %                                                padding 'my_padding'.
    %
    % Author(s): Aidan Gerkis
    % Date: 04-04-2024
    
    % Arrays of allowed input values
    in_names = ["im", "kernel", "padding"];
    in_types = ["numeric", "numeric", "char"];
    padding_allowed = ["symmetric", "replicate", "circular"];
    
    argin = cell(1, nargin); % Store input arguments for parsing

    switch nargin % Process inputs & assign default values if needed
        case 1 % Not enough input arguments, throw error
            kernel = fspecial('gaussian', [9, 9], 4);
            padding = 'circular';
            argin{1} = im;
        case 2 % Need to assign padding type
            padding = 'circular';
            argin{1} = im;
            argin{2} = kernel;
        case 3 % All arguments passed
            argin{1} = im;
            argin{2} = kernel;
            argin{3} = padding;
        otherwise
            error("Error in function call, too many input arguments");
    end

    for i=1:nargin % Check input types
        if ~isa(argin{i}, in_types(i))
            error("Unexpected type for parameter '" + in_names(i) + ...
                "'. Expected a " + in_types(i) + " but got a " + class(argin{i}) + ".");
        end
    end

   if ~ismember(padding, padding_allowed) % Validate input types
        error("Unknown padding type requested.");
   end

   b = imfilter(im, kernel, padding, 'conv'); % Filter image
end