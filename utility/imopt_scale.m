% imopt_scale.m
%
% Reads an input image, converts it to Black & White, and limits the pixel
% range to [0, 1]. Note that the image must be on your MATLAB path for this
% function to work. Implements code provided by Courtney Paquette in [1].
%
% Inputs:
%   im_name: The image name. [String]
%
% Outputs:
%   im: The image, converted to Black & White and scaled. [m x n Matrix]
%
% Author(s): Aidan Gerkis
% Date: 02-04-2024
%
% References:
%   [1]: C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" 
%        in MATH 564 - Honours Convex Optimization.

function im = imopt_scale(im_name)
    I = im2gray(imread(im_name)); % Load image as black & white
    
    % Scale image
    I = double(I(:, :, 1));
    mn = min(I(:));
    I = I - mn;
    mx = max(I(:));
    im = I/mx; % Scale image so each pixel is in [0, 1]
end