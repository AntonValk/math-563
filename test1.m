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

%im_clean = Primal_DouglasRachford_Splitting(b, kernel, 1, 0.5, 0.05, 250);
im_clean = PrimalDual_DouglasRachford_Splitting(b, kernel, 0.5, 0.25, 0.01, 250);
figure('Name', "Image after de-blurring");
imshow(im_clean, []);

%% Super Naive Parameter Sweep
%step_sizes = [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0];
% step_sizes = [0.5, 1.0, 1.5];
% reg_const = [0.1, 0.25, 0.5];
% gammas = [0.01, 0.025, 0.05];
% 
% for i=1:length(step_sizes)
%     for j=1:length(reg_const)
%         for k=1:length(gammas)
%             im_clean = PrimalDual_DouglasRachford_Splitting(b, kernel, step_sizes(i), reg_const(j), gammas(k), 50);
%             figure('Name', "Image after de-blurring "+ num2str(step_sizes(i)) + ": " + num2str(reg_const(j))  + ": " + num2str(gammas(k)));
%             imshow(im_clean, []);
%         end
%     end
% end
