function xSol = Primal_DouglasRachford_Splitting(b, t, rho)
[numRows, numCols]=size(b);
z1 = zeros(numRows, numCols);   %z1 is for boxProx
z21 = zeros(numRows, numCols);  % z21 is for l1/l2^2norm
z22 = zeros(numRows, numCols);  %z22 and z23 are for isoNorm
z23 = zeros(numRows, numCols);
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
t=2.0    %t = i.tprimaldr; % Need to change this for the various algorithms you are applying
applyMat = @(x) x + applyKTrans(applyK(x)) + applyDTrans(applyD(x));
eigValsMat = ones(numRows, numCols) + t*t*eigArry_KTrans.*eigArry_K + t*t*eigArry_D1Trans.*eigArry_D1...
    + t*t*eigArry_D2Trans.*eigArry_D2;

%R^(m x n) Computing (I + K^T*K + D^T*D)^(-1)*x
invertMatrix = @(x) ifft2(fft2(x)./eigValsMat);
for k = 1:500
    w=zeros(numRows,numCols,2);
    w(:,:,1)=z22;   % put z22 and z23 into a 3d matrix of size nxnx2
    w(:,:,2)=z23;
    x=boxProx(z1);
    y1=l1Prox(z21); %y1 is a matrix of size nxn from l1Prox
    y2=isoProx(w);
    aux = 2*x(:) - z1(:)+(applyKTrans(2*y1(:)-z21(:))+applyDTrans(2*y2(:)-[z22;z23])); %x(:) and z(:) is to vectorize nxn matrix to a n^2 by 1 vector
    u=invertMatrix(aux);
    v=[applyK(u(1:numRows,:));applyD1(u(numRows+1:2*numRows,:));applyD2(u(2*numRows+1:end,:))];
    z1=z1+rho*(u-x);
    z21=z21+rho*(v(1:numRows,:)-y1);
    z22=z22+rho*(v(numRows+1:2*numRows,:)-y2(1:numRows,:));
    z23=z23+rho*(v(2*numRows+1:end,:)-y2(2*numRows+1:end,:));
end
xSol=boxProx([z1;z21;z22;z23]);
end 