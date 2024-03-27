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
im_clean = PrimalDual_DouglasRachford_Splitting(b, kernel, 0.5, 0.25, 0.01, 50);
%im_clean = chambolle_pock(b, kernel, 1, 1, 0.05, 15);

figure('Name', "Image after de-blurring");
imshow(im_clean, []);

%% Super Naive Parameter Sweep
step_sizes = [0.01, 0.025, 0.04, 0.05, 0.1];
reg_const = [0.1, 0.25, 0.5, 1.5];
ss = [0.001, 0.025, 0.5, 1, 10]; %, 0.04, 0.05, 0.1];
ts = ss; 
gammas = [0.25, 0.5]; %[0.001, 0.0025, 0.005];

for i=1:length(reg_const) %length(ss)
    % Initialize plots for reg_const(i)
    figure('Name', "Image after de-blurring - Regularization: " + num2str(reg_const(i)));
    t = tiledlayout(length(step_sizes), length(gammas));

    for j=1:length(step_sizes) %length(ts)
        for k=1:length(gammas)
            % De-blur image
            im_clean = PrimalDual_DouglasRachford_Splitting(b, kernel, step_sizes(j), reg_const(i), gammas(k), 5);
            %im_clean = chambolle_pock(b, kernel, ss(i), ts(j), gammas(k), 25);
            %figure('Name', "Image after de-blurring "+ num2str(ss(i)) + ": " + num2str(ts(j))  + ": " + num2str(gammas(k)));
            
            nexttile; % Plot deblurred image
            imshow(im_clean, []);

            if j==1 % Add plot titles to top of layout
                title("\gamma = " + num2str(gammas(k)),'FontWeight','Normal');
            end

            if k==1 % Add plot titles to left side of layout
                ylabel("t = " + num2str(step_sizes(j)));
            end
        end
    end
end
