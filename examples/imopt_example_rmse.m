% imopt_example_rmse.m
%
% An example showing how to compute the error using RMSE.

clear; clc; close all;

I_original = imopt_scale('cameraman.jpg'); % Convert to B&W on the range [0, 1]
figure('Name','image before blurring')
imshow(I_original,[])

% Apply blur and noise
[I_blurred, kernel] = imopt_corrupt(I_original); % Apply blur & noise

figure('Name','image after blurring')
imshow(I_blurred,[]);

% Define custom parameters
p = struct();
p.metric = 'rmse'; % Use rmse to compute error
p.x_true = I_original; % Measure error w.r.t the original image
p.display = true; % Make plots

[I_deblurred, ~, D] = imopt(I_blurred, kernel, 'primal_dr', p); % Deblur image

figure('Name', "Image after de-blurring"); % Show deblurred image
imshow(I_deblurred, []);