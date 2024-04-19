# IMOPT
Welcome to IMOPT! IMOPT is a non-blind image deblurring package designed and implemented by Aidan Gerkis, April Niu, Antonios Valkanas, Cheng Shou, Linda Hu
 for the MATH 563 Final Project. This package implements four different deblurring algorithms, the Primal Douglas-Rachford algorithm, the Primal-Dual Douglas Rachford algorithm, the Alternating Direction Method of Multipliers algorithm, and the Chambolle-Pock algorithm. For those interested, more details on these algorithms and the theory behind them may be found in [1] & [2], but for now lets jump right into using IMOPT to deblur your images!

# Requirements
To use IMOPT we recommend running a version of MATLAB from r2022a onwards.

# Installation
To start download the IMOPT code. This can be done by cloning the repo or downloading the .zip containing the code (click the big green button in the top right). If you choose to download the .zip then place the zip in the final directory where you want to install IMOPT and unzip it. You should now have a folder in which you can see all of the project files. Open MATLAB and navigate to the directory where you placed the repo or unzipped the code too. Right-click on the file `imopt_install.m` and select "Run". This will set up your file path and run tests to ensure all the code works properly. A console output will inform you of the test status, if you see any errors in this log check your MATLAB version and contact us for support!

When you see the output
```
IMOPT test completed unsuccessfully. Passed 88 of 88 tests.
```
IMOPT has been successfully installed, and you can move on to deblurring your pictures!

# Usage
## Package Overview
The function `imopt` performs the heavy-lifting of deblurring the image, taking as an input the blurred image, the kernel of the blur, and some optional parameter. It then outputs the deblurred image, as well as some optional outputs with more information about the deblurring process. This function takes the following form:
```
[im_clean, error_final, D] = imopt(b, kernel, alg, p) 
```
The function supports four possible inputs, two of which are mandatory.
- `b`: The image to be deblurred. This image must be a 2D-matrix of black & white pixels, normalized to be in the range [0, 1]. This input is mandatory!
- `kernel`: The kernel used to blur the image. This should be formatted as a 2D-matrix and is mandatory!
- `alg`: The deblurring algorithm to use when deblurring the image, formatted as a character array. This input is optional, and if it is not specified then the default algorithm, Primal Douglas-Rachford, is used.
- `p`: A structure containing parameters for the deblurring algorithm.
 
The function also supports anywhere from 1 to 3 outputs, which are
 - `im_clean`: The deblurred image. This output is black & white, with pixel values in the range [0, 1].
 - `loss_final`: The final value of the loss function.
 - `D`: The complete output structure of the algorithm, containing detailed information on the optimization procedure.

To corrupt an image with blur and noise you can use the function `imopt_corrupt`, which applies a default blur and noise. For more control over the type of blur and noise the functions `imopt_blur` and `imopt_noisify` can be called directly.

## Quick Start
To get started with image deblurring first place your desired image anywhere on the MATLAB path. For these examples we will use "cameraman.jpg", a default image available in the \images folder. First call `imopt_scale` to ensure your image is black & white with pixels in the range [0, 1]
```
I_original = imopt_scale('cameraman.jpg', 256);
```
the second input ensures that the largest image dimension is 256 pixels, which ensures fast computations. You can then blur the image using 'imopt_corrupt'
```
[I_blurred, kernel] = imopt_corrupt(I_original);
```
where the first output is the blurred image and the second output is the kernel used to blur the image. We are then ready to call imopt!
```
I_deblurred = imopt(I_blurred, kernel)
```
To display the deblurred image we can use MATLAB's native 'imshow'
```
figure;
imshow(I_deblurred, []);
```
and voila, you have deblurred your first image! The image will also be saved to your active MATLAB directory.

## Advanced Parameter Options
The `imopt` function can be called with many different parameters, the first of which is the `alg` parameter, which specifies the deblurring algorithm to use. This function may take one of four values:
- `'primal_dr'`: Uses the Primal Douglas-Rachford algorithm to deblur the image.
- `'primaldual_dr'`: Uses the Primal-Dual Douglas-Rachford algorithm to deblur the image.
- `'admm'`: Uses the ADMM algorithm to deblur the image.
- `'chambolle_pock'`: Uses the Chambolle-Pock algorithm to deblur the image.

If the `alg` parameter is specified without specifying `p` then the default algorithm parameters for the chosen algorithm are used. However, the parameter structure `p` can be exploited for more control over the algorithm performance. This parameter structure accepts a wide selection of optional inputs.
- `x_init`: The initial guess to use when starting the algorithms. [An m x n matrix] [Default: All zeros]
- `x_true`: The true image. Mandatory for error calculations. [An m x n matrix] [Default: Empty]
- `t`: The step-size parameter for all algorithms. [A double] [Default: See report]
- `s`: The step-size parameter for Chambolle-Pock. [Double] [Default: See report]
- `g`: The denoising constant for all algorithms. [Double] [Default: See report]
- `rho`: The regularization constant for Primal Douglas-Rachford, Primal-Dual Douglas-Rachford, and ADMM. [Double] [Default: See report]
- `regularization`: The deblurring term to use. ['L1' or 'L2'] [Default: 'L1']
- `metric`: The error metric to use when evaluating error. ['rmse', 'psnr', or 'varinfo'] [Default: 'rmse']
- `L`: The number of bins to use when evaluating the variation of information. Mandatory if 'metric'=='varinfo'. [Double] [Default: N/A]
- `max_iters`: The maximum number of iterations to run the algorithm for. [Double] [Default: 500]
- `e_t`: The threshold error for the early stop criteria. [Double] [Default: 0.01]
- `tol`: The threshold loss for the early stop criteria. [Double] [Default: 0]
- `verbose`: A boolean indicating the level of console output. [Logical] [Default: false]
- `silent`: A boolean indicating if all console output should be stopped. [Logical] [Default: false]
- `display`: A boolean indicating if performance plots should be made. [Logical] [Default: true]
- `save_iters`: An integer indicating the level of saving of algorithm iterates. [0: No saving, 1: Full Saving, 2: Sparse Saving] [Default: 0]
- `ns`: The step size to use with sparse saving. [Integer] [Default:50]
- `im_name`: The name under which to save the deblurred image. [String] [Default: ""]
- `dir`: The directory under which to save the deblurred image. [String] [Default: current directory]

As an example, let's run the Chambolle-Pock algorithm using the `'L2'` deblurring term and custom parameter values. To compute the error we will use RMSE and provide the true image to the parameter structure. We'll also set an early stop threshold and turn on verbose mode. First, we need to make the parameter structure
```
p = struct();
p.regularization = 'L2';
p.t = 0.1;
p.s = 0.15;
p.g = 0.5;
p.tol = 0.5;
p.verbose = true;
p.x_true = I_original;
```
we can now call the 'imopt' function on the image we had previously blurred
```
[I_deblurred, loss_final] = imopt(I_blurred, kernel, 'chambolle_pock', p);
```

## Advanced Output Options
The `imopt` function supports three outputs, the third of which is the output structure from the respective algorithm calls, containing detailed information on the optimization process. The parameters in this structure, `D`, are summarized below
- `xf`: The final image after deblurring. [An m x n matrix]
- `t`: The time it took to run the optimization algorithm. [Double, Seconds]
- `k_end`: The number of iterations of the algorithm. [Integer]
- `e_end`: The error at the final iteration. [Double]
- `ek`: The error at each iteration. [1 x k_end Matrix]
- `fk`: The loss function value at each iteration. [1 x k_end Matrix]
- `xk`: The image at each iteration. Only saved if 'save_iters > 0'. [m x n x k_end Matrix]
- `inputs`: The parameter structure passed to the algorithm functions by 'imopt'. [Struct]

## Advanced Blurring Options
Custom blurs and noises can be applied using the `imopt_blur` & `imopt_noisify` functions. Lets use these features to apply a motion blur and salt & pepper noise with 50% density to the "cameraman.jpg" image. First we need to load and scale the image
```
I_original = imopt_scale("cameraman.jpg");
```
to apply a blur we first must create a kernel using the MATLAB `fspecial` function (https://www.mathworks.com/help/images/ref/fspecial.html). This kernel can then be used with `imopt_blur` to blur the image. The `imopt_blur` function supports three different types of boundary conditions, "circular", "replicate", and "symmetric". By default the "circular" boundary conditions are used, although in this example we apply the "replicate" conditions.
```
kernel = special('motion', 10, 15);
I_blurred = imopt_blur(I_original, kernel, 'replicate');
```
We can then apply noise to the image using `imopt_noisify`, which supports the same inputs as MATLAB's native `imnoise` function (https://www.mathworks.com/help/images/ref/imnoise.html)
```
I_noised = imopt_noisify(I_blurred, 'salt & pepper', 0.5);
```
we can now deblur the image using `imopt`
```
I_deblurred = imopt(I_noised, kernel);
```
Note that we recommend applying blur first and noise second, as this better represents the physical processes generating these phenomenon. Many different types of kernels and noise are supported, see the MATLAB documentation for more information.

## Advanced Display Options
The IMOPT package provides a function, `imopt_display`, with which convergence rate plots can be generated for a given IMOPT output structure `D`. The function takes three inputs
- `D`: The output structure from a call to IMOPT. [Struct]
- `plot`: The plot to generate. ['Error Evolution', 'Convergence', 'Loss Evolution', 'Loss Convergence', 'Image Iterates']
- 'n': The number of times to display the anmiation generated by 'Image Iterates'. [Double]
The plots generated by 'Error Evolution' and 'Convergence' plot the error versus iterations and log of error versus iterations, respectively, and thus require that the true image be provided when calling `imopt`. Examples of these output plots are shown below

<img src = "https://github.com/AntonValk/math-563/assets/44247293/fadbe185-da10-4e43-90cc-d6974dbcf3e4" width="350" height="350">
<img src = "https://github.com/AntonValk/math-563/assets/44247293/c66ec812-68ee-438e-9237-680cd6916a12" width="350" height="350">

The plots generated by 'Loss Evolution' and 'Loss Convergence plot the loss versus iterations and log of loss versus iterations, respectively. Example outputs are shown below

<img src = "https://github.com/AntonValk/math-563/assets/44247293/59449191-d4d8-434b-b58c-414111cf2bed" width="350" height="350">
<img src = "https://github.com/AntonValk/math-563/assets/44247293/593e9f7f-1ed7-49df-a903-7e28d781390a" width="350" height="350">

The option 'Image Iterates' plots the image guess at each iteration in an animation and requires `save_iters > 0`. To generate this animation we first must set up our parameter structure appropriately
```
p = struct();
p.save_iters = 2; % Enables sparse saving, better for memory performance
p.ns = 20; % Save every 20th iterate
```
we can then call `imopt` on a blurred image
```
I_original = imopt_scale('cameraman.jpg');
[I_blurred, kernel] = imopt_corrupt(I_original);
[I_deblurred, loss_final, D] = imopt(I_blurred, kernel, 'primal_dr', p);
```
we can then call `imopt_display` on the output structure `D` to create the animation
```
imopt_display(D, 'Image Iterates', 1);
```
and finally we can save the image iterates as .pngs.
```
imopt_save_iterates(D);
```
Note that calling this function in this way, without specifying the directory to save the images, saves the images in the current directory.

# Conclusion
Hopefully you are now ready to begin your image deblurring journey! Thank you for choosing IMOPT and please contact us with any questions or concerns!

Best,

Aidan Gerkis, April Niu, Antonios Valkanas, Cheng Shou, Linda Hu

# References

[1] C. Paquette, "MATH 463/563 - Convex Optimization, Project Description" in MATH 564 - Honours Convex Optimization.

[2] Daniel O’Connor and Lieven Vandenberghe. Primal Dual Decomposition by Operator Splitting and Applications to Image Deblurring. Siam Journal Imaging Sciences, 7(3), 2014.
