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

%im_clean = primal_douglasrachford_splitting(b, kernel, 1, 0.5, 1.5, 50); %0.05, 250);
%im_clean = primaldual_douglasrachford_splitting(b, kernel, 0.005, 0.005, 0.1, 150); % <- In general smaller parameters seem to work better here, sweep those
im_clean = admm(b, kernel, 0.005, 0.005, 0.1, 150);
%im_clean = chambolle_pock(b, kernel, 1e-12, 1e-12, 0.5, 250);

figure('Name', "Image after de-blurring");
imshow(im_clean, []);

%% Super Naive Parameter Sweep
% step_sizes = [0.0001, 0.0005, 0.001, 0.005, 0.01];
% reg_const = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1];
% ss = [0.000001, 0.00001, 0.00005, 0.0001, 0.0005, 0.001];
% ts = ss; 
% gammas = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1];
% 
% for i=1:length(gammas) %length(reg_const) %length(ss)
%     % Initialize plots for reg_const(i)
%     %figure('Name', "Image after de-blurring - Regularization: " + num2str(reg_const(i)));
%     figure('Name', "Image after de-blurring - Gamma: " + num2str(gammas(i)));
%     %t = tiledlayout(length(step_sizes), length(gammas));
%     t = tiledlayout(length(ss), length(ts));
% 
%     for j=1:length(ss)%length(step_sizes) %length(ts)
%         for k=1:length(ts)
%             % De-blur image
%             %im_clean = PrimalDual_DouglasRachford_Splitting(b, kernel, step_sizes(j), reg_const(i), gammas(k), 50);
%             im_clean = chambolle_pock(b, kernel, ss(j), ts(k), gammas(i), 25);
%             %figure('Name', "Image after de-blurring "+ num2str(ss(i)) + ": " + num2str(ts(j))  + ": " + num2str(gammas(k)));
%             
%             nexttile; % Plot deblurred image
%             imshow(im_clean, []);
% 
%             if j==1 % Add plot titles to top of layout
%                 %title("\gamma = " + num2str(gammas(k)),'FontWeight','Normal');
%                 title("s = " + num2str(ss(k)),'FontWeight','Normal');
%             end
% 
%             if k==1 % Add plot titles to left side of layout
%                 %ylabel("t = " + num2str(step_sizes(j)));
%                 ylabel("t = " + num2str(ts(j)));
%             end
%         end
%     end
% end