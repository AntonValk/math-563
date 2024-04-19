% Test image
I_original = imopt_scale('cameraman.jpg');

% Specify blurred image and kernel
[b, kernel] = imopt_corrupt(I_original); % ------------------------SPECIFY BLURRED IMAGE AND KERNEL HERE
x = zeros(256, 256); % <---------------------------------------------------ADJUST STARTING POINT HERE

% Create parameter structure
params = struct();
params.x_init = x;
params.verbose = true; % Enable print outputs
params.display = true;

I_deblurred_pdr = imopt(b, kernel, 'primal_dr', params); % Deblur image with Primal DR

I_deblurred_pddr = imopt(b, kernel, 'primaldual_dr', params); % Deblur image with Primal DR

I_deblurred_admm = imopt(b, kernel, 'admm', params); % Deblur image with Primal DR

I_deblurred_cp = imopt(b, kernel, 'chambolle_pock', params); % Deblur image with Primal DR

% Show images
figure('Name', 'Deblurred: Primal DR');
imshow(I_deblurred_pdr, []);

figure('Name', 'Deblurred: Primal-Dual DR');
imshow(I_deblurred_pddr, []);

figure('Name', 'Deblurred: ADMM');
imshow(I_deblurred_admm, []);

figure('Name', 'Deblurred: Chambolle-Pock');
imshow(I_deblurred_cp, []);