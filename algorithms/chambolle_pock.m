function  xSol=chambolle_pock(b, kernel, t, s, g, k_max)
[numRows, numCols]=size(b);

%computes the numRow x numCol matrix of the eigenvalues for K and D1 and
%D2; Here D1 = I oplus D1 in the paper and D2 = D1 oplus I.
eigArry_K = eigValsForPeriodicConvOp(kernel, numRows, numCols);
eigArry_D1 = eigValsForPeriodicConvOp([-1,1]', numRows, numCols);
eigArry_D2 = eigValsForPeriodicConvOp([-1,1], numRows, numCols);

%computes numRow x numCol matrix of the eigenvalues for K^T and D1^T and
%D2^T;
eigArry_KTrans = conj(eigArry_K);
eigArry_D1Trans = conj(eigArry_D1);
eigArry_D2Trans = conj(eigArry_D2);

%Functions which compute Kx, D1x, D2x, Dxt, K^Tx, D1^Tx, D2^Tx, and D^Ty.
%Note for all the x functions, the input x is in R^(m x n) and outputs into
%R^(m x n) except for D which outputs into 2 concat. R^(m x n) matrices;
%For D^Ty, y is two m x n matrices concatanated and outputs into R^(m x n)
applyK = @(x) applyPeriodicConv2D(x, eigArry_K);
applyD1 = @(x) applyPeriodicConv2D(x, eigArry_D1);
applyD2 = @(x) applyPeriodicConv2D(x, eigArry_D2);

applyKTrans = @(x) applyPeriodicConv2D(x, eigArry_KTrans);
applyD1Trans = @(x) applyPeriodicConv2D(x, eigArry_D1Trans);
applyD2Trans = @(x) applyPeriodicConv2D(x, eigArry_D2Trans);

applyD = @(x) cat(3, applyD1(x), applyD2(x));

applyDTrans = @(y) applyD1Trans(y(:,:,1)) + applyD2Trans(y(:, :, 2));

% Function which computes the (I + K^TK + D^TD)x where x in R^(m x n)
% matrix and the eigenvalues of I + t*t*K^TK + t*t*D^TD; here t is the
% stepsizes
applyMat = @(x) x + applyKTrans(applyK(x)) + applyDTrans(applyD(x));
eigValsMat = ones(numRows, numCols) + t*t*eigArry_KTrans.*eigArry_K + t*t*eigArry_D1Trans.*eigArry_D1...
    + t*t*eigArry_D2Trans.*eigArry_D2;

%R^(m x n) Computing (I + K^T*K + D^T*D)^(-1)*x
invertMatrix = @(x) ifft2(fft2(x)./eigValsMat);

% Initialize
xk = 0.1 + zeros(numRows, numCols); % n x n matrix
zk = zeros(numRows, numCols); % n x n matrix
yk = [applyK(xk); applyD1(xk); applyD2(xk)]; % 3n x n matrix

for k = 1:k_max
    xk_old = xk; % Save previous xk

    % Compute Prox Ops
    wk = yk + s*[applyK(zk); applyD1(zk); applyD2(zk)]; % Input to prox of g
    wk = reshape(wk, [numRows, numCols, 3]);
    
    yk = reshape(yk, [numRows, numCols, 3]);
    vk = xk - t*(applyKTrans(yk(:, :, 1)) + applyD1Trans(yk(:, :, 2)) + applyD2Trans(yk(:, :, 3))); % Input to prox of f

    % Prox of sg*
    yk1 = wk(:, :, 1) - l1Prox(wk(:, :, 1) - b, s) - b; % Part one of prox of sg*
    yk2 = wk(:, :, 2:3) - isoProx(wk(:, :, 2:3), s*g); % Part two of prox of sg*
    yk = [yk1; yk2(:, :, 1); yk2(:, :, 2)]; % Compile components of prox of sg*

    % Prox of tf
    xk = boxProx(vk);

    % Update zk
    zk = 2*xk - xk_old;
end
xSol = xk;
end