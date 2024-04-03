%   hist2.m 
%   Calculates the joint histogram of two images or signals
%
%   Usage:
%   n=hist2(A,B,L) is the joint histogram of matrices A and B, using L
%   bins for each matrix.
%
%   Inputs:
%   A: A matrix in R^(n x n) that represents the ground truth image.
%   B: A matrix in R^(n x n) that represent the reconstructed image.
%   L: Number of bins of discrete distribution.
%
%   Author: Antonios Valkanas
%   Date: 02-04-2024
% 
%   Adapted from: https://www.mathworks.com/matlabcentral/fileexchange/20688-kullback-leibler-divergence

function n = hist2(A, B, L)
ma=min(A(:)); 
MA=max(A(:)); 
mb=min(B(:)); 
MB=max(B(:));
% For sensorimotor variables, in [-pi,pi] 
% ma=-pi; 
% MA=pi; 
% mb=-pi; 
% MB=pi;
% Scale and round to fit in {0,...,L-1} 
A=round((A-ma)*(L-1)/(MA-ma+eps)); 
B=round((B-mb)*(L-1)/(MB-mb+eps)); 
n=zeros(L); 
x=0:L-1; 
for i=0:L-1 
    n(i+1,:) = histc(B(A==i),x,1); 
end
end
