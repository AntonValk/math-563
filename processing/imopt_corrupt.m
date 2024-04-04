% imopt_corrupt.m
%
% Applies default corruption to an image.
%
% Inputs:
%   im: The image to corrupt. [m x n Matrix]
%
% Outputs:
%   c: The corrupted image. [m x n Matrix]
%   k: The kernel used to blur the image. [k x k Matrix]
%
% Author(s): Aidan Gerkis
% Date: 04-04-2024

function [c, k] = imopt_corrupt(im)
    if ~isa(im, "numeric") % Validate input type
        error("Error in function call, expected a numeric but got a " + class(im));
    end
    
    k = fspecial('gaussian', [9, 9], 4); % Define default kernel
    b = imopt_blur(im, k); % Blur
    c = imopt_noisify(b); % Add noise
end