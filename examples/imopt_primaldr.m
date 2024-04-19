% test_imopt

clear; clc; close all;

I = imopt_scale('cameraman.jpg'); % Convert to B&W on the range [0, 1]
figure('Name','image before blurring')
imshow(I,[])

% Apply blur and noise
[b, kernel] = imopt_corrupt(I); 
figure('Name','image after blurring')
imshow(b,[]);

im_in = cat(3, b, I); % Concatenate with true image

% Create parameter structure
params = struct();
params.regularization = 'L1'; % Select between L1 & L2 regularization
params.t = 0.35; % Define custom step size
params.rho = 0.5; % Define custom regularization constant
params.gamma = 0.25; % Define custom gamma
params.e_t = 0.01; % Set early stop based on error threshold
params.tol = 1; % Set early stop based on loss
params.max_iter = 50; % Sets maximum number of iterations
params.verbose = true; % Print additional information about package status
params.display = true; % Display images summarizing algorithm performance
params.save_iters = 1; % Save all iterates at each algorithm step

[im_clean, ~, D] = imopt(im_in, kernel, 'primal_dr', params); % Deblur image

figure('Name', "Image after de-blurring"); % Show deblurred image
imshow(im_clean, []);