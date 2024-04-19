% imopt_example_primaldualdr.m
%
% An example showing how to use Primal-Dual Douglas-Rachford to deblur an image.

clear; clc; close all;

I_original = imopt_scale('cameraman.jpg'); % Convert to B&W on the range [0, 1]
figure('Name','image before blurring')
imshow(I_original,[])

% Apply blur and noise
[I_blurred, kernel] = imopt_corrupt(I_original); 
figure('Name','image after blurring')
imshow(I_blurred,[]);

% Create parameter structure
p = struct();
p.regularization = 'L1'; % Select between L1 & L2 regularization
p.t = 1.6; % Define custom step size
p.rho = 1.5; % Define custom regularization constant
p.gamma = 0.001; % Define custom gamma
p.e_t = 0.01; % Set early stop based on error threshold
p.tol = 1; % Set early stop based on loss
p.max_iter = 50; % Sets maximum number of iterations
p.verbose = true; % Print additional information about package status
p.display = true; % Display images summarizing algorithm performance
p.save_iters = 1; % Save all iterates at each algorithm step
p.x_true = I_original; % True image

I_deblurred = imopt(I_blurred, kernel, 'primaldual_dr', p); % Deblur image

figure('Name', "Image after de-blurring"); % Show deblurred image
imshow(I_deblurred, []);