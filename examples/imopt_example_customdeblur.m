% imopt_example_customdeblur
%
% An example showing the user how to create an image with a custom blur
% applied

clear; clc; close all;

I_original = imopt_scale('cameraman.jpg'); % Convert to B&W on the range [0, 1]
figure('Name','image before blurring')
imshow(I_original,[])

% Apply blur and noise
kernel = fspecial('motion', 10, 15); % Create kernel
I_blurred = imopt_blur(I_original, kernel, 'symmetric'); % Apply blur
I_noised = imopt_noisify(I_original, 'salt & pepper', 0.5); % Apply noise

figure('Name','image after blurring')
imshow(I_noised,[]);

I_deblurred = imopt(I_noised, kernel); % Deblur image

figure('Name', "Image after de-blurring"); % Show deblurred image
imshow(I_deblurred, []);