I = imread('cameraman.jpg');
I = rgb2gray(I);   %black and white
I = double(I(:, :, 1));
mn = min(I(:));
I = I - mn;
mx = max(I(:));
I = I/mx;
figure('Name','image before blurring')
imshow(I,[])
kernel = fspecial('gaussian', [15,15], 5);
b = imfilter(I,kernel);
noiseDensity=0.5;
b = imnoise(b,'salt & pepper',noiseDensity); 
figure('Name','image after blurring')
imshow(b,[]);
Primal_DouglasRachford_Splitting(b, 1, 0.5);
