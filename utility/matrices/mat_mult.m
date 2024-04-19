function Ax = mat_mult(x, A, K, t)
    % mat_mult.m
    %
    % A function wrapper that calls the matrix multiplication code written by
    % Courtney Paquette and distributed as part of the assignment.
    %
    % Inputs:
    %   x: The matrix to multiply. [matrix or tensor, depending on desired operation]
    %   A: The multiplication to perform. [char array]
    %   K: The image kernel. [matrix]
    %   t: A constant to scale the matrix multiplication, only needed for
    %      computing the inverse. [Double]
    %
    % Outputs:
    %   Ax: The results of the matrix multiplication. [matrix or tensor, depending on desired operation]
    %
    % Author(s): Aidan Gerkis
    % Date: 01-04-2024
    
    % Preliminaries - Computing eigenvalues
    [numRows, numCols]=size(x);

    %computes the numRow x numCol matrix of the eigenvalues for K and D1 and
    %D2; Here D1 = I oplus D1 in the paper and D2 = D1 oplus I.
    eigArry_K = eigValsForPeriodicConvOp(K, numRows, numCols);
    eigArry_D1 = eigValsForPeriodicConvOp([-1,1]', numRows, numCols);
    eigArry_D2 = eigValsForPeriodicConvOp([-1,1], numRows, numCols);
    
    %computes numRow x numCol matrix of the eigenvalues for K^T and D1^T and
    %D2^T;
    eigArry_KTrans = conj(eigArry_K);
    eigArry_D1Trans = conj(eigArry_D1);
    eigArry_D2Trans = conj(eigArry_D2);

    switch A % Perform multiplication
        case 'K' % Input is R^(m x n), Output is R^(m x n)
            Ax = applyPeriodicConv2D(x, eigArry_K);
        case 'D1' % Input is R^(m x n), Output is R^(m x n)
            Ax = applyPeriodicConv2D(x, eigArry_D1);
        case 'D2' % Input is R^(m x n), Output is R^(m x n)
            Ax = applyPeriodicConv2D(x, eigArry_D2);
        case 'KT' % Input is R^(m x n), Output is R^(m x n)
             Ax = applyPeriodicConv2D(x, eigArry_KTrans);
        case 'D1T' % Input is R^(m x n), Output is R^(m x n)
            Ax = applyPeriodicConv2D(x, eigArry_D1Trans);
        case 'D2T' % Input is R^(m x n), Output is R^(m x n)
            Ax = applyPeriodicConv2D(x, eigArry_D2Trans);           
        case 'D' % Input is R^(m x n), Output is 2 concatenated R^(m x n) matrices
            Ax = cat(3, mat_mult(x, 'D1', K), mat_mult(x, 'D2', K));
        case 'DT' % Input is 2 concatenated R^(m x n) matrices, Output is R^(m x n)
            Ax = mat_mult(x(:, :, 1), 'D1T', K) + mat_mult(x(:, :, 2), 'D2T', K);
        case 'ATA' % Computes (I + K^TK + D^TD)x, Input is R^(m x n), Output is R^(m x n)
            Ax = x + mat_mult(mat_mult(x, 'K', K), 'KT', K) + mat_mult(mat_mult(x, 'D', K), 'DT', K);
        case 'inv' % Computes (I + t*t*K^T*K + t*t*D^T*D)^(-1)*x, Input is R^(m x n), Output is R^(m x n)
            eigValsMat = ones(numRows, numCols) + t*t*eigArry_KTrans.*eigArry_K + t*t*eigArry_D1Trans.*eigArry_D1...
                + t*t*eigArry_D2Trans.*eigArry_D2;
            Ax = ifft2(fft2(x)./eigValsMat);
        otherwise
            error("Unknown matrix multiplication requested.");
    end
end