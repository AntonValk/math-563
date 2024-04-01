% test_imopt

clear; clc; close all;

I = imread('cameraman.jpg');
I = rgb2gray(I);   %black and white
I = double(I(:, :, 1));
mn = min(I(:));
I = I - mn;
mx = max(I(:));
I = I/mx;
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
params.max_iter = 50; % Maximum number of iterations
params.verbose = true; % Print additional information about function run
params.display = true; % Display images
params.save_iters = true; % Save iterates

%[im_clean, ~, D] = imopt(b, kernel, 'primal_dr', params);
[im_clean, ~, D] = imopt(b, kernel, 'primaldual_dr', params);

figure('Name', "Image after de-blurring");
imshow(im_clean, []);