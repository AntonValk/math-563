function xSol = Primal_DouglasRachford_Splitting(b, kernel, t, rho)
[numRows, numCols]=size(b);
% z1 = zeros(numRows, numCols);   %z1 is for boxProx
% z21 = zeros(numRows, numCols);  % z21 is for l1/l2^2norm
% z22 = zeros(numRows, numCols);  %z22 and z23 are for isoNorm
% z23 = zeros(numRows, numCols);

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
%z1 = b;
z1 = zeros(numRows, numCols);
z2 = [applyK(z1); applyD1(z1); applyD2(z1)];

for k = 1:5 %%TODO: Better stopping condition lol <- The paper suggests a stopping condition on the norm of the resolvents acting on x (pg5)
    % Split z2 up
    z21 = z2(1:numRows, :);
    z22 = cat(3, z2((numRows + 1):2*numRows, :), z2((2*numRows + 1):end, :)); % Change z2 from 2D (2n*n)->3D (n*n*2)
    
    % Compute prox ops
    x = boxProx(z1);
    y1 = l1Prox(z21 - b, t) + b; % From Equation 8 in the reference
    y2 = isoProx(z22, t);

    % Compute resolvent of B
    appK = applyKTrans(2*y1 - z21); % Compute the K component of the multiplication
    appD = applyDTrans(2*y2 - z22); % Compute the D component of the multiplication (applyDTrans uses the concatenated form of the arrays)
    u = invertMatrix(2*x - z1 + appK + appD); % Compute the resolvent
    
    % Compute Updates
    v = [applyK(u); applyD1(u); applyD2(u)];
    y = [y1; y2(:, :, 1); y2(:, :, 2)]; % Need to build y as a 2d matrix
    
    z1 = z1 + rho*(u - x);
    z2 = z2 + rho*(v - y);
end
xSol=boxProx(z1);
end 