function  PrimalDual_DouglasRachford_Splitting(b, t, rho)
[numRows, numCols]=size(b);
p = zeros(numRows, numCols);
q1 = zeros(numRows, numCols);  % q1 is for l1/l2^2norm
q2 = zeros(numRows, numCols);  % q2 and q3 are for isoNorm
q3 = zeros(numRows, numCols);  
t=1.0  %t = i.tprimaldr; % Need to change this for the various algorithms you are applying
%computes the numRow x numCol matrix of the eigenvalues for K and D1 and
%D2; Here D1 = I oplus D1 in the paper and D2 = D1 oplus I.
kernel = fspecial('gaussian', [15,15], 5);
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
t=1.0  %t = i.tprimaldr; % Need to change this for the various algorithms you are applying
applyMat = @(x) x + applyKTrans(applyK(x)) + applyDTrans(applyD(x));
eigValsMat = ones(numRows, numCols) + t*t*eigArry_KTrans.*eigArry_K + t*t*eigArry_D1Trans.*eigArry_D1...
    + t*t*eigArry_D2Trans.*eigArry_D2;

%R^(m x n) Computing (I + K^T*K + D^T*D)^(-1)*x
invertMatrix = @(x) ifft2(fft2(x)./eigValsMat);
for k = 1:500
    x = boxProx(p,t);   %x is n by n matrix
    w=zeros(numRows,numCols,2);
    w(:,:,1)=q2;   % put z22 and z23 into a 3d matrix of size nxnx2
    w(:,:,2)=q3;
    y1 = l1Prox(q1,t);
    y2 = isoProx(w,t);
    y = [y1;y2(:,:,1);y2(:,:,2)];   %y is 3n by n matrix
    z = [q1;q2:q3]-y;   %z is 3n by n matrix
    inverse = invertMatrix(eye(numRows^2));
    iTA = [eye(numRows^2);t*applyK(eye(numRows^2));t*applyD1(eye(numRows^2));t*applyD2(eye(numRows^2))];
    iTATrans = [eye(numRows^2),-t*applyKTrans(eye(numRows^2)),-t*applyD1Trans(eye(numRows^2));-t*applyD2Trans(eye(numRows^2))];
    bigProd = iTA.*inverse.*iTATrans;
    blockWithId = [zeros(2*numRows^2),zeros(2*numRows^2);zeros(2*numRows^2),eye(2*numRows^2)];
    bigSum = blockWithId+bigProd;
    aux = [2*x-p;2*z-[q1;q2:q3]];
    res = bigSum .* aux(:);
    w = reshape(res(1:numRos^2,:),numRows,numCols);
    v = reshap(res(numRows^2+1:end,:),3*numRows, numCols);
    p = p+rho*(w-x);
    q = q+rho*(v-z);
end
xSol = boxProx(p,t);
end