% imopt_example_deblur
%
% A script showcasing the basics of how to blur and deblur an image with
% imopt.
%
% Author: Aidan Gerkis
% Date: 18-04-2024

clear; clc; close all;

I = imopt_scale('cameraman.jpg'); % Convert to B&W on the range [0, 1]
figure('Name','image before blurring')
imshow(I,[])

% Apply blur and noise
[b, kernel] = imopt_corrupt(I); 
figure('Name','image after blurring')
imshow(b,[]);

im_in = cat(3, b, I); % Concatenate with true image

[im_clean, ~, D] = imopt(im_in, kernel); % Deblur image

figure('Name', "Image after de-blurring"); % Plot deblurred image
imshow(im_clean, []);