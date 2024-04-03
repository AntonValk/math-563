% mutual_info.m
%
% Inputs:
%   A: A matrix in R^(n x n) that represents the ground truth image.
%   B: A matrix in R^(n x n) that represent the reconstructed image.
%   L: Number of bins of histogram.
%
% Outputs:
%   I: Mutual Information (MI) (scalar)
%
% Usage:
%   I = mutual_info(A, B, L) outputs the mutual information between 
%   the reconstructed image and the ground truth
%
%   Author: Antonios Valkanas
%   Date: 02-04-2024
%
%   MUTUALINFO Calculates the mutual information between two images.
%
%   I = MUTUALINFO(A, B)        Mutual information of images A and B of
%   identical size, using L bins for histograms.
%
%   Adapted from <a href="matlab:web('https://www.mathworks.com/matlabcentral/fileexchange/13289-fast-mutual-information-of-two-images-or-signals')">MI</a> by Jose Delpiano.
%
%   References:
%   [1] T. M. Cover and J. A. Thomas, "Entropy, Relative Entropy, and
%       Mutual Information," in Elements of Information Theory, 2nd ed.
%       Hoboken, NJ: Wiley-Interscience, 2006, pp. 13–55.


function I=mutual_info(A,B,varargin) 
if nargin>=3
    L=varargin{1};
else
    L=256;
end
A=double(A); 
B=double(B); 
     
na = hist(A(:),L); 
na = na/sum(na);
nb = hist(B(:),L); 
nb = nb/sum(nb);
n2 = hist2(A,B,L); 
n2 = n2/sum(n2(:));
I=sum(minf(n2,na'*nb)); 
% -----------------------
function y=minf(pab,papb)
I=find(papb(:)>1e-12 & pab(:)>1e-12); % function support 
y=pab(I).*log2(pab(I)./papb(I));
