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
im_clean = PrimalDual_DouglasRachford_Splitting(b, kernel, 1, 0.05, 0.25, 15);

figure('Name', 'Image after de-blurring');
imshow(im_clean, []);