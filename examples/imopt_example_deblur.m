% imopt_example_deblur
%
% A script showcasing the basics of how to blur and deblur an image with
% imopt.
%
% Author: Aidan Gerkis
% Date: 18-04-2024

clear; clc; close all;

I_original = imopt_scale('cameraman.jpg'); % Convert to B&W on the range [0, 1]
figure('Name','image before blurring')
imshow(I_original,[])

% Apply blur and noise
[I_blurred, kernel] = imopt_corrupt(I_original); 
figure('Name','image after blurring')
imshow(I_blurred,[]);

I_deblurred = imopt(I_blurred, kernel); % Deblur image

figure('Name', "Image after de-blurring"); % Plot deblurred image
imshow(I_deblurred, []);