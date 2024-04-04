% test_imopt

clear; clc; close all;

I = imopt_scale('cameraman.jpg'); % Convert to B&W on the range [0, 1]
figure('Name','image before blurring')
imshow(I,[])

% Apply blur and noise
kernel = fspecial('gaussian', [9,9], 4);
%kernel = fspecial('motion', 50, 45);
b = imfilter(I,kernel);
noiseDensity=0.1;
b = imnoise(b, 'salt & pepper', noiseDensity); 
figure('Name','image after blurring')
imshow(b,[]);

im_in = cat(3, b, I);
%im_in = padarray(im_in, [9,9], 'both'); % Add zero padding

% Create parameter structure
params = struct();
params.regularization = 'L1';
%params.s = 0.35;
%params.t = 0.35; % Step size
%params.gamma = 0.25; % Gamma
%params.rho = 0.5; % Regularization constant
params.e_t = 0.01; % Error threshold
params.tol = 1; % Loss threshold
params.max_iter = 150; % Maximum number of iterations
params.verbose = true; % Print additional information about function run
params.display = true; % Display images
params.save_iters = false; % Save iterates

[im_clean, ~, D] = imopt(im_in, kernel, 'primal_dr', params);
%[im_clean, ~, D] = imopt(im_in, kernel, 'primaldual_dr', params);
%[im_clean, ~, D] = imopt(im_in, kernel, 'admm', params);
%[im_clean, ~, D] = imopt(im_in, kernel, 'chambolle_pock', params);

figure('Name', "Image after de-blurring");
imshow(im_clean, []);