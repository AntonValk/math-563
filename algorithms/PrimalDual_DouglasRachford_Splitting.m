function  xSol=PrimalDual_DouglasRachford_Splitting(b, kernel, t, rho, g, k_max)
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
pk = zeros(numRows, numCols); % n x n matrix
qk= zeros(3*numRows, numCols);  % 3n x n matrix

for k = 1:k_max
    % Extract elements
    q1 = qk(1:numRows, :);
    q2 = cat(3, qk((numRows + 1):2*numRows, :), qk((2*numRows + 1):end, :)); % Change q from 2D (2n*n)->3D (n*n*2)

    % Compute Prox Ops
    xk = boxProx(pk);   %x is n by n matrix
    z1 = l1Prox(q1 - b, t) + b; % Equation 8 in reference
    z2 = isoProx(q2, t*g);
    zk = [q1; q2(:, :, 1); q2(:, :, 2)] - [z1; z2(:, :, 1); z2(:, :, 2)];

    % Compute Resolvent of B (pg 7 in reference) <- TODO: Double check
    inverse = invertMatrix(ones(numRows, 1)); % inverse(I + t^2A'A)
    iTA = [eye(numRows); t*applyK(eye(numRows)); t*applyD1(eye(numRows)); t*applyD2(eye(numRows))]; % [I; tA]
    iTATrans = [eye(numRows); -t*applyK(eye(numRows)); -t*applyD1(eye(numRows)); -t*applyD2(eye(numRows))]'; % [I, -tA']
    bigProd = iTA*inverse*iTATrans; % Combine above

    blockWithId = [zeros(3*numRows, 4*numRows); zeros(numRows, 3*numRows), eye(numRows)]; % <- Also not confident that the dimensions on this are right!
    %blockWithId = [zeros(2*numRows), zeros(2*numRows); zeros(2*numRows), eye(2*numRows)];

    res_b = blockWithId + bigProd; % Compile resolvent

    res_b_z = res_b*[2*xk - pk; 2*zk - qk]; % Compute the vector satisfing 0 in res_b
    
    wk = res_b_z(1:numRows, :); % In n x n
    vk = res_b_z((numRows + 1):end, :); % In 3n x n

    % Perform update
    pk = pk + rho*(wk - xk); % n x n matrix
    qk = qk + rho*(vk - zk); % 3n x n matrix
end
xSol = boxProx(pk);
end
