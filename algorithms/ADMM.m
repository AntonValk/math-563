function xSol = ADMM(b, kernel, t, rho, g, k_max)
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

% Initialize xk, uk, wk, yk, zk
%z1 = b;
xk = zeros(numRows, numCols); %xk, uk, wk are n by n matrices
uk = zeros(numRows, numCols);
wk = zeros(numRows, numCols);
y1 = zeros(numRows, numCols);%yk,zk are 3n by n matrices
y2 = zeros(numRows, numCols);
y3 = zeros(numRows, numCols);
z1 = zeros(numRows, numCols);
z2 = zeros(numRows, numCols);
z3 = zeros(numRows, numCols);


for k = 1:k_max %%TODO: Better stopping condition lol <- The paper suggests a stopping condition on the norm of the resolvents acting on x (pg5)
    atz = applyKTrans(z1) + applyD1Trans(z2) + applyD2Trans(z3); %(A^t z)
    aty = applyKTrans(y1) + applyD1Trans(y2) + applyD2Trans(y3); %(A^t y)
    xk = invertMatrix(uk + aty - (1/t)*(wk + atz)); %
    % z2 = mat_split(z2k, 3);
    %z21 = z2k(1:numRows, :);
    %z22 = cat(3, z2k((numRows + 1):2*numRows, :), z2k((2*numRows + 1):end, :)); % Change z2 from 2D (2n*n)->3D (n*n*2)
    
    % Compute prox ops
    uk = boxProx(rho*xk + (1-rho)*uk + wk/t);
    w = mat_split([applyK(xk);applyD1(xk); applyD2(xk)],3);
    yk = mat_split([y1; y2; y3], 3);
    zk = mat_split([z1; z2; z3], 3);
    % unsure about the +/-b and */t for the two prox operators.
    y1 = l1Prox(t*(rho* w(:,:,1) + (1-rho)*y1+z1/t)+b,t)/t-b; % From Equation 8 in the reference 
    y_aux = isoProx(rho*w(:,:,2:3) + (1-rho)*yk(:,:,2:3) + zk(:,:,2:3)/t, g*t)/t;
    %yk = [y1; y_aux(:,:,1); y_aux(:,:,2)];
    
    % Compute Updates
    wk = wk + t*(xk-uk);
    z1 = z1 + t * (applyK(xk) - y1);
    z2 = z2 + t * (applyD1(xk) - y_aux(:,:,1));
    z3 = z3 + t * (applyD2(xk) - y_aux(:,:,2));
end
atz = applyKTrans(z1) + applyD1Trans(z2) + applyD2Trans(z3); %(A^t z)
aty = applyKTrans(y1) + applyD1Trans(y2) + applyD2Trans(y3); %(A^t y)
xSol= invertMatrix(uk + aty - (1/t)*(wk + atz));
end 

