function imopt_install()
    % imopt_install.m
    %
    % Checks version and adds current path to MATLAB path, effectively
    % installing IMOPT.
    %
    % Author: Aidan Gerkis
    % Date: 19-04-2024

    % Check version
    out_of_date = ~isMATLABReleaseOlderThan("R2022b"); % Recommended to run on r2022a or later
    cur_version = matlabRelease;
    cur_release = cur_version.Release;

    if out_of_date % Print warning if version is out of date
        fprintf("Your MATLAB release is currently %s, we recommend to run IMOPT on MATLAB R2022a or later.\n", cur_release);
    end

    % Add path
    cur_dir = pwd; % Get current directory
    path = genpath(cur_dir); % Find all folders in current directory
    addpath(path); % Update path

    % Test install
    imopt_test;
end