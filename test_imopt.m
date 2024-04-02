% test_imopt

clear; clc; close all;

I = imopt_scale('cameraman.jpg'); % Convert to B&W on the range [0, 1]
figure('Name','image before blurring')
imshow(I,[])
kernel = fspecial('gaussian', [15,15], 5);
b = imfilter(I,kernel);
noiseDensity=0.1;
b = imnoise(b, 'salt & pepper', noiseDensity); 
figure('Name','image after blurring')
imshow(b,[]);

% Create parameter structure
params = struct();
params.regularization = 'L1';
params.e_t = 0.001; % Error threshold
params.max_iter = 500; % Maximum number of iterations
params.verbose = true; % Print additional information about function run
params.display = true; % Display images
params.save_iters = true; % Save iterates

im_in = cat(3, b, I);
[im_clean, ~, D] = imopt(im_in, kernel, 'primal_dr', params);
%[im_clean, ~, D] = imopt(b, kernel, 'primaldual_dr', params);
%[im_clean, ~, D] = imopt(b, kernel, 'admm', params);
%[im_clean, ~, D] = imopt(b, kernel, 'chambolle_pock', params);

figure('Name', "Image after de-blurring");
imshow(im_clean, []);