function [true, inputs] = imopt_get_true()
    % imopt_get_true.m
    %
    % Returns the expected correct values to use when testing the IMOPT
    % package. True results are stored in the test_results.m &
    % test_results.mat files.
    %
    % Author: Aidan Gerkis
    % Date: 19-04-2024
    
    % Load true results
    test_results;
    % Initiate Outputs
    true = struct();
    inputs = struct();

    % Files to Check
    true.imopt_files = ["admm", "chambolle_pock", "primal_douglasrachford_splitting",...
                        "primaldual_douglasrachford_splitting", "parallel_sweep", ...
                        "parallel_sweep_benchmark", "conv_rate", "deblur_movie", ...
                        "error_evo", "imopt_display", "loss_conv_rate", "loss_evo", ...
                        "imopt_example_deblur", "imopt_primaldr", "cameraman.jpg", ...
                        "manWithHat.tiff", "mcgill.jpg", "imopt_inds", "imopt_iso", ...
                        "imopt_l1", "imopt_l2", "imopt_loss", "imopt_psnr", ...
                        "imopt_rmse", "imopt_var_info", "imopt_blur", "imopt_corrupt", ...
                        "imopt_noisify", "imopt_scale", "boxProx", "isoProx", "l1Prox", ...
                        "l2squaredProx", "prox_ln", "prox_tg", "benchmark_parallel", ...
                        "early_stop", "hist2", "imopt_get_default", "imopt_save_iterates", ...
                        "imopt_test", "memory_check", "mutual_info", "set_plotting_parameters", ...
                        "eigValsForPeriodicConvOp", "applyPeriodicConv2D", "applyXTrans", ...
                        "mat_mult", "mat_split", "imopt"];

    % Results of matrix tests
    inputs.x = x;
    inputs.w = w;
    inputs.A = [x; x];
        inputs.numRows = nr;
    inputs.numCols = nc;
    inputs.filter = kernel;
    inputs.kernelsize = length(kernel(1, :));

    inputs.k = 2;
    inputs.t = 1;
    inputs.mult = ["K", "D1", "D2", "KT", "D1T", "D2T", "D", "DT", "ATA", "inv"];

    true.eVal = eVals; % Eigenvalue test results

    true.axconv = axconv; % Convolution test results
    true.atxconv = atxconv; % Transpose test results

    true.split = matrix_split; % mat_split test results
    true.ax = struct(); % mat_mult test results
    true.ax.K = Kx;
    true.ax.D1 = D1x;
    true.ax.D2 = D2x;
    true.ax.KT = KTx;
    true.ax.D1T = D1Tx;
    true.ax.D2T = D2Tx;
    true.ax.D = Dx;
    true.ax.DT = DTx;
    true.ax.ATA = ATAx;
    true.ax.inv = invx;

    % Results of prox tests
    inputs.b = b;
    
    true.box = boxp;
    true.l1 = l1_prox;
    true.l2 = l2_prox;
    true.iso = iso_prox;
    true.ln1 = ln1;
    true.ln2 = ln2;
    true.tg_ln1 = tg_ln1;
    true.tg_ln2 = tg_ln2;

    % Results of image processing tests
    inputs.im_name = im_name;

    true.blur = blur;
    true.noise = noise;
    true.scale = scale;
    
    % Results of metric tests
    inputs.nbins = nbins;

    true.hist2 = binned;
    true.info = mut_info;
    true.vi = var_info;
    true.rmse = rootmean_err;
    true.psnr = psn_ratio;
    true.ind = ind;
    true.iso_norm = iso_norm;
    true.l1_norm = l1_norm;
    true.l2_norm = l2_norm;
    true.loss_test = loss;

    % Results of algorithm tests
    inputs.x_init = x_init;
    inputs.pdrt = pdrt;
    inputs.pdrg = pdrg;
    inputs.pdrrho = pdrrho;
    inputs.fpdr = fpdr;
    inputs.pddrt = pddrt;
    inputs.pddrg = pddrg;
    inputs.pddrrho = pddrrho;
    inputs.fpddr = fpddr;
    inputs.admmt = admmt;
    inputs.admmg = admmg;
    inputs.admmrho = admmrho;
    inputs.fadmm = fadmm;
    inputs.cpt = cpt;
    inputs.cps = cps;
    inputs.cpg = cpg;
    inputs.fcp = fcp;
    inputs.k_max = k_max;
    inputs.e_t = e_t;

    true.pdr = pdr_res;
    true.pddr = pddr_res;
    true.admm = admm_res;
    true.cp = cp_res;

    % Results of IMOPT Package test
    true.im_clean = im_clean;