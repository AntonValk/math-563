% generate_test_results.m
%
% Generate the 'true' results on which to test the IMOPT install.
%
% Author: Aidan Gerkis
% Date: 19-04-2024

clc; clear;

x = imopt_scale('cameraman.jpg', 256);
b = imopt_corrupt(x);

w = x;
[nr, nc] = size(x);

% Matrix tests
kernel = fspecial('gaussian', [9, 9], 4);
eVals = eigValsForPeriodicConvOp(kernel, nr, nc);
axconv = applyPeriodicConv2D(x, eVals);
atxconv = applyXTrans(w, x, length(kernel(1, :)));
matrix_split = mat_split([x; x], 2);
Kx = mat_mult(x, 'K', kernel, 1);
D1x = mat_mult(x, 'D1', kernel, 1);
D2x = mat_mult(x, 'D2', kernel, 1);
KTx = mat_mult(x, 'KT', kernel, 1);
D1Tx = mat_mult(x, 'D1T', kernel, 1);
D2Tx = mat_mult(x, 'D2T', kernel, 1);
Dx = mat_mult(x, 'D', kernel, 1);
DTx = mat_mult(cat(3, x, w), 'DT', kernel, 1);
ATAx = mat_mult(x, 'ATA', kernel, 1);
invx = mat_mult(x, 'inv', kernel, 1);

% Prox tests
boxp = boxProx(x, 1);
l1_prox = l1Prox(x, 1);
l2_prox = l2squaredProx(x, 1);
iso_prox = isoProx(cat(3, x, w), 1);
ln1 = prox_ln(x, 1, 'L1', b);
ln2 = prox_ln(x, 1, 'L2', b);
tg_ln1 = prox_tg(cat(3, x, x, x), 1, 1, 'L1', b);
tg_ln2 = prox_tg(cat(3, x, x, x), 1, 1, 'L2', b);

% Image processing tests
im_name = 'cameraman.jpg';
rng(1);
blur = imopt_blur(x, kernel);
rng(1);
noise = imopt_noisify(x);
scale = imopt_scale(im_name);

% Metric tests
nbins = 50;

binned = hist2(x, b, nbins);
mut_info = mutual_info(x, b, nbins);
var_info = imopt_var_info(x, b, nbins);
rootmean_err = imopt_rmse(x, b);
psn_ratio = imopt_psnr(x, b);
ind = imopt_inds(x);
iso_norm = imopt_iso(cat(3, x, w));
l1_norm = imopt_l1(x);
l2_norm = imopt_l2(x);
loss = imopt_loss(x, b, 1, 'L1', kernel);

% Algorithm tests
x_init = zeros(256, 256);
functions = struct();
functions.err_eval = @(x_hat) imopt_rmse(x, x_hat);
functions.prox_l = @(x, t) prox_ln(x, t, 'L1', b);
functions.early_stop = @(cur, prev) early_stop(cur, prev, 0.1);
pdrt = 0.1;
pdrg = 0.01;
pdrrho = 1.05;
pddrt = 1.75;
pddrg = 0.05;
pddrrho = 1.75;
admmt = 0.7407;
admmg = 1e-15;
admmrho = 1.4815;
cpt = 0.1;
cps = 0.4072;
cpg = 1e-15;
k_max = 50;
e_t = 0.001;
fpdr = functions;
fpdr.f_loss = @(x)imopt_loss(x, b, pdrg, 'L1', kernel);
fpddr = functions;
fpddr.f_loss = @(x)imopt_loss(x, b, pdrg, 'L1', kernel);
fadmm = functions;
fadmm.f_loss = @(x)imopt_loss(x, b, pdrg, 'L1', kernel);
fcp = functions;
fcp.f_loss = @(x)imopt_loss(x, b, pdrg, 'L1', kernel);

pdr_res = primal_douglasrachford_splitting(b, kernel, x_init, fpdr, pdrt, pdrg, pdrrho, k_max, e_t, 0, 0, 0);
pddr_res = primaldual_douglasrachford_splitting(b, kernel, x_init, fpddr, pddrt, pddrg, pddrrho, k_max, e_t, 0, 0, 0);
admm_res = admm(b, kernel, x_init, fadmm, admmt, admmg, admmrho, k_max, e_t, 0, 0, 0);
cp_res = chambolle_pock(b, kernel, x_init, fcp, cpt, cps, cpg, k_max, e_t, 0, 0, 0);

% IMOPT Package tests
im_clean = imopt(b, kernel);