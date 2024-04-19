% imopt_example_imageiterates_animation.m
%
% An example showing how to save image iterates, create an animation, and
% save iterates as pngs.

clear; clc; close all;

I = imopt_scale('cameraman.jpg'); % Convert to B&W on the range [0, 1]
figure('Name','image before blurring')
imshow(I,[])

% Apply blur and noise
[I_blurred, kernel] = imopt_corrupt(I); % Apply blur & noise

figure('Name','image after blurring')
imshow(I_blurred,[]);

% Define custom parameters
p = struct();
p.save_iters = 2; % Sparse save of iterates
p.ns = 20; % Save every 20 iterates
p.display = false; % Turn displays off

[I_deblurred, ~, D] = imopt(I_blurred, kernel, 'primal_dr', p); % Deblur image

figure('Name', "Image after de-blurring"); % Show deblurred image
imshow(I_deblurred, []);

% Create an animation and play it five times
imopt_display(D, 'Image Iterates', 5);

% Save image iterates in outputs folder
cur_folder = cd;
folder = fullfile(cur_folder, "..");
dir = folder + "\outputs";

imopt_save_iterates(D, dir);