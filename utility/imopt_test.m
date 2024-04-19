function imopt_test()
    % imopt_test.m
    %
    % Tests important functions in IMOPT, ensuring functions execute successfully
    % and checking results to ensure they are correct.
    %
    % Author: Aidan Gerkis
    % Date: 19-04-2024

    [t, in] = imopt_get_true; % Get true test results
    
    % Status variables
    passed_tot = 0;
    total_tot = 88;

    % Initial messge
    fprintf("Checking install of IMOPT. If tests are skipped then fix any errors encountered and rerun imopt_test.\n");

    % Check for all files in path
    fprintf("Checking MATLAB path...");
    passed = 0; % Indicates number of passed tests
    file_failed = 0; % Flag to indicate if file checking failed
    idx = 1;
    
    errors = struct(); % Store error messages
    warnings = struct(); % Store warnings
    
    for i=1:length(t.imopt_files) % Iterate over all files
        file_exists = exist(t.imopt_files(i), 'file'); % Check existence
    
        if file_exists ~= 2 % 2 indicates that the file exists
            errors.(char(idx + 64)) = "Could not find " + t.imopt_files(i) + ...
                " on MATLAB path.";
            file_failed = 1;
        else % Update number of passed tests
            passed = passed + 1;
        end
    end
    
    passed_tot = passed_tot + passed;

    % Display status
    fprintf("passed %d of %d\n", passed, length(t.imopt_files));
    display_errors_warnings(errors, warnings);
    
    % Check that matrix code works
    passed = 0; % Indicates number of passed tests
    mult_failed = 0; % Flag to indicate if matrix multiplication checks failed
    total = 14; % Total number of tests
    idx = 1;
    
    errors = struct(); % Store error messages
    warnings = struct(); % Store warnings
    
    fprintf("Testing matrix multiplications...")
    
    if ~file_failed % Only run if previous tests passed
        try % Check eigenvalue calculations
            eVal_test = eigValsForPeriodicConvOp(in.filter, in.numRows, in.numCols);
        
            if ~isequal(eVal_test, t.eVal) % Display error if result doesn't match ground truth
                errors.(char(idx + 64)) = "Error in eigValsForPeridodicConvOp.m, incorrect output. Check ifft2 functionality.";
                mult_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            mult_failed = 1;
            idx = idx + 1;
        end
    else % Set warning if tests were skipped
        warnings.(char(idx + 64)) = "Warning: skipped tests of mat_split.m, mat_mult.m.";
        idx = idx + 1;
    end
    
    
    if ~mult_failed && ~ file_failed % Only run if previous tests passed
        try % Check matrix multiplication
            ax_test = applyPeriodicConv2D(in.x, eVal_test);
    
            if ~isequal(ax_test, t.axconv) % Display error if result doesn't match ground truth
                errors.(char(idx + 64)) = "Error in applyPeriodicConv2D.m, incorrect output. Check ifft2 functionality.";
                mult_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            mult_failed = 1;
            idx = idx + 1;
        end
    
    
        try % Check transpose matrix multiplication
            atx_test = applyXTrans(in.w, in.x, in.kernelsize);
    
            if ~isequal(atx_test, t.atxconv) % Display error if result doesn't match ground truth
                errors.(char(idx + 64)) = "Error in applyXTrans.m, incorrect output. Check conv2 functionality.";
                mult_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            mult_failed = 1;
            idx = idx + 1;
        end
    else % Set warning if tests were skipped
        warnings.(char(idx + 64)) = "Warning: skipped tests of applyPeriodicConv2D.m, applyXTrans.m.";
        idx = idx + 1;
    end
    
    if ~mult_failed && ~ file_failed % If previous tests passed
        try % Test matrix splitting
            split_test = mat_split(in.A, in.k);

            if ~isequal(split_test, t.split) % Display error if result doesn't match ground truth
                errors.(char(idx + 64)) = "Error in mat_split.m, incorrect output. Check cat functionality.";
                mult_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            mult_failed = 1;
            idx = idx + 1;
        end
    
        fields = fieldnames(t.ax);
        for i=1:10 % Test matrix multiplications
            try % Test matrix
                if strcmp(fields{i}, 'DT') % Two possible inputs to mat_mult
                    ax_test = mat_mult(cat(3, in.x, in.w), in.mult(i), in.filter, in.t);
                else
                    ax_test = mat_mult(in.x, in.mult(i), in.filter, in.t);
                end
    
                if ~isequal(ax_test, t.ax.(fields{i})) % Check result of multiplication
                    errors.(char(idx + 64)) = "Error in mat_mult.m, incorrect output for case " + in.mult(i) + ".";
                    mult_failed = 1;
                    idx = idx + 1;
                else
                    passed = passed + 1;
                end
            catch err % Save error and set flag to skip next tests if failed
                errors.(char(idx + 64)) = getReport(err);
                mult_failed = 1;
                idx = idx + 1;
            end
        end
    else % Set warning if tests were skipped
        warnings.(char(idx + 64)) = "Warning: skipped tests of mat_split.m, mat_mult.m.";
    end
    
    passed_tot = passed_tot + passed;

    % Display test results
    fprintf("passed %d of %d.\n", passed, total);
    display_errors_warnings(errors, warnings);
    
    % Test prox operators
    passed = 0; % Indicates number of passed tests
    prox_failed = 0; % Flag to indicate if matrix multiplication checks failed
    total = 6; % Total number of tests
    idx = 1;
    
    errors = struct(); % Store error messages
    warnings = struct(); % Store warnings
    
    fprintf("Testing proximal operators...")
    
    if ~file_failed % Only run if all files exist
        try % Test box prox
            box_test = boxProx(in.x, 1);

            if ~isequal(box_test, t.box) % Check result
                errors.(char(idx + 64)) = "Error in boxProx.m, incorrect output.";
                prox_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            prox_failed = 1;
            idx = idx + 1;
        end

        try % Test L1 Prox
            l1_test = l1Prox(in.x, 1);

            if ~isequal(l1_test, t.l1) % Check result
                errors.(char(idx + 64)) = "Error in l1Prox.m, incorrect output.";
                prox_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            prox_failed = 1;
            idx = idx + 1;
        end

        try % Test L2 Prox
            l2_test = l2squaredProx(in.x, 1);

            if ~isequal(l2_test, t.l2) % Check result
                errors.(char(idx + 64)) = "Error in l2squaredProx.m, incorrect output.";
                prox_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            prox_failed = 1;
            idx = idx + 1;
        end

        try % Test ISO Prox
            iso_test = isoProx(cat(3, in.x, in.w), 1);

            if ~isequal(iso_test, t.iso) % Check result
                errors.(char(idx + 64)) = "Error in isoProx.m, incorrect output.";
                prox_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            prox_failed = 1;
            idx = idx + 1;
        end
        
        if ~prox_failed % Only proceed if all previous proximal operator tests were passed
            try % Test ln Prox
                l1_test = prox_ln(in.x, 1, 'L1', in.b);
                l2_test = prox_ln(in.x, 1, 'L2', in.b);

                if ~isequal(l1_test, t.ln1) || ~isequal(l2_test, t.ln2) % Check result
                    errors.(char(idx + 64)) = "Error in prox_ln.m, incorrect output.";
                    prox_failed = 1;
                    idx = idx + 1;
                else
                    passed = passed + 1;
                end
            catch err % Save error and set flag to skip next tests if failed
                errors.(char(idx + 64)) = getReport(err);
                prox_failed = 1;
                idx = idx + 1;
            end

            try % Test prox_tg
                l1_test = prox_tg(cat(3, in.x, in.x, in.x), 1, 1, 'L1', in.b);
                l2_test = prox_tg(cat(3, in.x, in.x, in.x), 1, 1, 'L2', in.b);

                if ~isequal(l1_test, t.tg_ln1) || ~isequal(l2_test, t.tg_ln2) % Check result
                    errors.(char(idx + 64)) = "Error in prox_tg.m, incorrect output.";
                    prox_failed = 1;
                else
                    passed = passed + 1;
                end
            catch err % Save error and set flag to skip next tests if failed
                errors.(char(idx + 64)) = getReport(err);
                prox_failed = 1;
            end
        else % Set warning if tests were skipped
            warnings.(char(idx + 64)) = "Warning: skipped tests of all prox_ln.m, prox_tg.m.";
        end
    else % Set warning if tests were skipped
        warnings.(char(idx + 64)) = "Warning: skipped tests of all proximal operators.";
    end
    
    passed_tot = passed_tot + passed;

    % Display test results
    fprintf("passed %d of %d.\n", passed, total);
    display_errors_warnings(errors, warnings);
    
    % Test image processing
    passed = 0; % Indicates number of passed tests
    process_failed = 0; % Flag to indicate if matrix multiplication checks failed
    total = 3; % Total number of tests
    idx = 1;
    
    errors = struct(); % Store error messages
    warnings = struct(); % Store warnings
    
    fprintf("Testing image processing...")
    
    if ~file_failed % Only run if all files exist
        try % Test blurring
            rng(1); % Force seed
            blur_test = imopt_blur(in.x, in.filter);

            if ~isequal(blur_test, t.blur) % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in imopt_blur.m. Check installation of image processing toolbox.";
                process_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            process_failed = 1;
            idx = idx + 1;
        end
        
        try % Test noising
            rng(1); % Force seed
            noise_test = imopt_noisify(in.x);

            if ~isequal(noise_test, t.noise) % Check result
                errors.(char(idx + 64)) = "Warning, incorrect output in imopt_noisify.m. Check installation of image processing toolbox.";
                process_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            process_failed = 1;
            idx = idx + 1;
        end

        try % Test scaling
            scale_test = imopt_scale(in.im_name);

            if ~isequal(scale_test, t.scale) % Check result
                errors.(char(idx + 64)) = "Warning, incorrect output in imopt_scale.m. Check installation of image processing toolbox.";
                process_failed = 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            process_failed = 1;
        end
    else % Set warning if tests were skipped
        warnings.(char(idx + 64)) = "Warning: skipped tests of image processing functions.";
    end
    
    passed_tot = passed_tot + passed;

    % Display test results
    fprintf("passed %d of %d.\n", passed, total);
    display_errors_warnings(errors, warnings);

    % Test metrics
    passed = 0; % Indicates number of passed tests
    metrics_failed = 0; % Flag to indicate if matrix multiplication checks failed
    total = 10; % Total number of tests
    idx = 1;
    
    errors = struct(); % Store error messages
    warnings = struct(); % Store warnings
    
    fprintf("Testing metric calculations...")
    
    if ~file_failed % Only run if all files exist
        try % Test hist2
            bin_test = hist2(in.x, in.b, in.nbins);

            if ~isequal(bin_test, t.hist2) % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in hist2.m. Check MATLAB version.";
                metrics_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            metrics_failed = 1;
            idx = idx + 1;
        end

        if ~metrics_failed % Don't test mutual info if hist2 failed
            try % Test mutual info
                info_test = mutual_info(in.x, in.b, in.nbins);
    
                if info_test ~= t.info % Check result
                    warnings.(char(idx + 64)) = "Warning, incorrect output in mutual_info.m.";
                    metrics_failed = 1;
                    idx = idx + 1;
                else
                    passed = passed + 1;
                end
            catch err % Save error and set flag to skip next tests if failed
                errors.(char(idx + 64)) = getReport(err);
                metrics_failed = 1;
                idx = idx + 1;
            end
        else % Set warning if tests were skipped
            warnings.(char(idx + 64)) = "Warning: skipped test of mutual_info.m.";
            idx = idx + 1;
        end

        if ~metrics_failed % Don't test mutual info if hist2 or mutual_info failed
            try % Test variation of information
                vi_test = imopt_var_info(in.x, in.b, in.nbins);
    
                if vi_test ~= t.vi % Check result
                    warnings.(char(idx + 64)) = "Warning, incorrect output in imopt_var_info.m.";
                    metrics_failed = 1;
                    idx = idx + 1;
                else
                    passed = passed + 1;
                end
            catch err % Save error and set flag to skip next tests if failed
                errors.(char(idx + 64)) = getReport(err);
                metrics_failed = 1;
                idx = idx + 1;
            end
        else % Set warning if tests were skipped
            warnings.(char(idx + 64)) = "Warning: skipped test of imopt_var_info.m.";
            idx = idx + 1;
        end
        
        try % Test rmse
            rmse_test = imopt_rmse(in.x, in.b);

            if rmse_test ~= t.rmse % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in imopt_rmse.m.";
                metrics_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            metrics_failed = 1;
            idx = idx + 1;
        end
        
        try % Test psnr
            psnr_test = imopt_psnr(in.x, in.b);

            if psnr_test ~= t.psnr % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in imopt_psnr.m.";
                metrics_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            metrics_failed = 1;
            idx = idx + 1;
        end
        
        test_loss = 1; % Flag to indicate if imopt_loss should be tested

        try % Test indicator
            ind_test = imopt_inds(in.x);

            if ind_test ~= t.ind % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in imopt_ind.m.";
                metrics_failed = 1;
                test_loss = 0;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            metrics_failed = 1;
            test_loss = 0;
            idx = idx + 1;
        end

        try % Test iso
            iso_test = imopt_iso(cat(3, in.x, in.w));

            if iso_test ~= t.iso_norm % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in imopt_iso.m.";
                metrics_failed = 1;
                test_loss = 0;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            metrics_failed = 1;
            test_loss = 0;
            idx = idx + 1;
        end

        try % Test l1 norm
            l1_test = imopt_l1(in.x);

            if l1_test ~= t.l1_norm % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in imopt_l1.m.";
                metrics_failed = 1;
                test_loss = 0;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            metrics_failed = 1;
            test_loss = 0;
            idx = idx + 1;
        end

        try % Test l2  norm
            l2_test = imopt_l2(in.x);

            if l2_test ~= t.l2_norm % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in imopt_l2.m.";
                metrics_failed = 1;
                test_loss = 0;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            metrics_failed = 1;
            test_loss = 0;
            idx = idx + 1;
        end

        if test_loss
            try % Test loss function
                loss_test = imopt_loss(in.x, in.b, 1, 'L1', in.filter);

                if loss_test ~= t.loss_test % Check result
                    warnings.(char(idx + 64)) = "Warning, incorrect output in imopt_loss.m.";
                    metrics_failed = 1;
                    idx = idx + 1;
                else
                    passed = passed + 1;
                end
            catch err % Save error and set flag to skip next tests if failed
                errors.(char(idx + 64)) = getReport(err);
                metrics_failed = 1;
                idx = idx + 1;
            end

        else % Set warning if tests were skipped
            warnings.(char(idx + 64)) = "Warning: skipped test of imopt_loss.m.";
            metrics_failed = 1;
            idx = idx + 1;
        end
    else % Set warning if tests were skipped
        warnings.(char(idx + 64)) = "Warning: skipped tests of metrics.";
        metrics_failed = 1;
    end
    
    passed_tot = passed_tot + passed;

    % Display test results
    fprintf("passed %d of %d.\n", passed, total);
    display_errors_warnings(errors, warnings);

    % Test algorithms
    passed = 0; % Indicates number of passed tests
    algs_failed = 0; % Flag to indicate if matrix multiplication checks failed
    total = 4; % Total number of tests
    idx = 1;
    
    errors = struct(); % Store error messages
    warnings = struct(); % Store warnings
    
    fprintf("Testing algorithms...")
    if ~metrics_failed && ~file_failed && ~prox_failed && ~mult_failed
        try % Test primal douglas-rachford
            pdr_test = primal_douglasrachford_splitting(in.b, in.filter, in.x_init, in.fpdr, in.pdrt, in.pdrg, in.pdrrho, in.k_max, in.e_t, 0, 0, 0);

            if ~isequal(pdr_test.xf, t.pdr.xf) % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in primal_douglasrachford_splitting.m.";
                algs_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            algs_failed = 1;
            idx = idx + 1;
        end

        try % Test primal-dual douglas-rachford
            pddr_test = primaldual_douglasrachford_splitting(in.b, in.filter, in.x_init, in.fpddr, in.pddrt, in.pddrg, in.pddrrho, in.k_max, in.e_t, 0, 0, 0);
        
            if ~isequal(pddr_test.xf, t.pddr.xf) % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in primaldual_douglasrachford_splitting.m.";
                algs_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            algs_failed = 1;
            idx = idx + 1;
        end

        try % Test admm
            admm_test = admm(in.b, in.filter, in.x_init, in.fadmm, in.admmt, in.admmg, in.admmrho, in.k_max, in.e_t, 0, 0, 0);

            if ~isequal(admm_test.xf, t.admm.xf) % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in admm.m.";
                algs_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            algs_failed = 1;
            idx = idx + 1;
        end

        try % Test chambolle-pock
            cp_test = chambolle_pock(in.b, in.filter, in.x_init, in.fcp, in.cpt, in.cps, in.cpg, in.k_max, in.e_t, 0, 0, 0);

            if ~isequal(cp_test.xf, t.cp.xf) % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in chambolle_pock.m.";
                algs_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            algs_failed = 1;
            idx = idx + 1;
        end
    else % Set warning if tests were skipped
        warnings.(char(idx + 64)) = "Warning: skipped tests of algorithms.";
        algs_failed = 1;
        idx = idx + 1;
    end
    
    passed_tot = passed_tot + passed;

    % Display test results
    fprintf("passed %d of %d.\n", passed, total);
    display_errors_warnings(errors, warnings);
    

    % Test imopt
    passed = 0; % Indicates number of passed tests
    imopt_failed = 0; % Flag to indicate success of imopt tests
    total = 1; % Total number of tests
    idx = 1;
    
    errors = struct(); % Store error messages
    warnings = struct(); % Store warnings
    
    fprintf("Testing IMOPT package...")
    if ~algs_failed && ~process_failed
        try % Test imopt
            params = struct();
            params.display = false;
            params.silent = true;
            im_test = imopt(in.b, in.filter, 'primal_dr', params);

            if ~isequal(im_test, t.im_clean) % Check result
                warnings.(char(idx + 64)) = "Warning, incorrect output in imopt.m.";
                imopt_failed = 1;
                idx = idx + 1;
            else
                passed = passed + 1;
            end
        catch err % Save error and set flag to skip next tests if failed
            errors.(char(idx + 64)) = getReport(err);
            imopt_failed = 1;
            idx = idx + 1;
        end
    else % Set warning if tests were skipped
        warnings.(char(idx + 64)) = "Warning: skipped test of IMOPT package.";
        imopt_failed = 1;
        idx = idx + 1;
    end
    
    passed_tot = passed_tot + passed;

    % Display test results
    fprintf("passed %d of %d.\n", passed, total);
    display_errors_warnings(errors, warnings);
    
    % Display final status
    success = ~metrics_failed && ~process_failed && ~file_failed && ~prox_failed && ~mult_failed && ~algs_failed && ~imopt_failed;
    if success
        status = "sucessfully";
    else
        status = "unsuccessfully";
    end

    fprintf("IMOPT test completed %s. Passed %d of %d tests.\n", status, passed_tot, total_tot);
end
