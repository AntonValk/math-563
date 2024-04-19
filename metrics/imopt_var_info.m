function I=imopt_var_info(A, B, L) 
    % imopt_var_info.m
    %
    % Inputs:
    %   A: A matrix in R^(n x n) that represents the ground truth image.
    %   B: A matrix in R^(n x n) that represent the reconstructed image.
    %   L: Number of bins of discrete distribution. Recommended value for a
    %      256x256 image is ~50.
    %
    % Outputs:
    %   I: Variation of information metric (VI) (scalar)
    %
    % Usage:
    %   I = imopt_var_info(A, B) outputs the mutual information between 
    %   the reconstructed image and the ground truth
    %
    % Author: Antonios Valkanas
    % Date: 02-04-2024
    %
    %   References:
    %   [1] T. M. Cover and J. A. Thomas, "Entropy, Relative Entropy, and
    %       Mutual Information," in Elements of Information Theory, 2nd ed.
    %       Hoboken, NJ: Wiley-Interscience, 2006.
    %   [2] https://en.wikipedia.org/wiki/Variation_of_information
    
    I = mutual_info(A, A, L) + mutual_info(B, B, L) + - 2 * mutual_info(A, B, L);
end