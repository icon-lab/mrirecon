%% Multilinear Coil Compression (MLCC) Demo
% This is a demo that shows an example implementation of MLCC algorithm
% over bSSFP data with multiple phase-cycled acquisitions and coils. It is
% based on Biyik et. al, "Reconstruction by Calibration over Tensors for
% Multi-Coil Multi-Acquisition Balanced SSFP Imaging". MLCC is a technique
% that enables the use of cs-MRI algorithms when disjoint sampling patterns
% are employed for bSSFP with multiple acquisitions. MLCC is based on
% multilinear singular value decomposition. In this demo, MLCC is compared
% to GCC (Zhang et. al, "Coil Compression for Accelerated Imaging with
% Cartesian Sampling") and ReCat algorithm is used for the reconstruction.


clearvars;
close all;

%% Parameters
N = 4; % number of phase-cycled acquisitions
R = 8; % undersampling factor (4,8,12,16 are present in the default package)
lambda = 0.018; % l_2 regularization term weight
beta = 0.05; % Tikhonov regularization parameter
kernelSize = [11,11]; % interpolation kernel size
calibSizePct = 0.13; % calibration region, percentage of the k-space size in
					 % both directions
nIter = 20; % number of LSQR iterations
wavWeight = 0; % Wavelet threshold term for l1 regularization (0 to disable)
p_acq = 4; % p value for p-norm combination over phase-cycles
p_coils = 2; % p value for p-norm combination over coils
Dp = 5; % number of virtual coils to be generated with GCC and MLCC
slwin = 1; % sliding window length for GCC (5 in ESPIRiT library GCC demo) 

%% Necessary libraries and folders
% SPIRiT V0.3 is required
if exist('CC.m', 'file')==0
	addpath(genpath('ESPIRiT'));
	if exist('CC.m', 'file')==0
		warning('ReCat utilizes SPIRiT V0.3 library. Please download from http://people.eecs.berkeley.edu/~mlustig/Software.html and copy the "ESPIRiT" directory to the same directory as this demo file. ');
		return;
	end
end
% TensorLab 3.0 is required
if exist('mlsvd.m', 'file')==0
	addpath(genpath('TensorLab'));
	if exist('mlsvd.m', 'file')==0
		warning('MLCC utilizes TensorLab 3.0 package. Please download from http://www.tensorlab.net/ and copy the "TensorLab" directory to the same directory as this demo file. ');
		return;
	end
end
% ReCat_code folder includes ReCat functions
if exist('recat_optimize.m', 'file')==0
	addpath('ReCat');
end
% util folder includes some utility functions
if exist('normalize.m', 'file')==0
	addpath('util');
end

%% DATA
% load sample bSSFP data
load('data/invivo_32coil.mat');
raw_data = double(raw_data); % LSQR implementation requires double type

%% Reference Image
% generate via p-norm combination over N=8, fully-sampled acquisitions
images = ifft2c(raw_data);
originalImage = sos(sos(images,4,p_coils),3,p_acq);

%% Undersampling Mask
% load pre-generated mask
load(['masks/mask_' num2str(R) 'x.mat']);

%% Prepare Data
% take N-many acquisitions
imageFFT = raw_data(:,:,1:N,:);
mask = mask(:,:,1:N);
% extract other parameters
D = size(images,4);
sx = size(images,1);
sy = size(images,2);
calibSize = round([calibSizePct*sx, calibSizePct*sy]);
% undersample
imageFFT = imageFFT .* repmat(mask, [1, 1, 1, D]);

%% Coil Compression
%% GCC
% GCC is applied separately to each phase-cycle, due to disjoint sampling
disp('GCC + ReCat is being performed');
imageFFT_gcc = zeros(sx, sy, N, Dp);
for acq_no=1:N
    acq = squeeze(imageFFT(:,:,acq_no,:));
    imageFFT_gcc(:,:,acq_no,:) = gcc_compress(acq, slwin, Dp); 
end

% ReCat over GCC compressed coils
% interpolation kernel estimation
disp('Estimating interpolation kernel - this may take a few minutes');
% crop the calibration region
kCalib = crop(imageFFT_gcc, [calibSize,N,Dp]);
% estimate the kernel over the tensor
kernel = zeros([kernelSize, N, Dp, N, Dp]);
YtY = data2YtY(kCalib, kernelSize);
for n=1:N
	for d=1:Dp
		kernel(:,:,:,:,n,d) = calibrate(YtY, kernelSize, N, Dp, n, d, beta);
	end
end
% LSQR reconstruction
disp('Performing LSQR reconstruction');
% perform the optimization
res = recat_optimize(imageFFT_gcc, kernel, nIter, lambda, wavWeight);
% transform to image domain
res_ims = ifft2c( res);
% P-norm combination
result_gcc = sos(sos(res_ims,4,p_coils),3,p_acq);


%% MLCC
% MLCC can be applied over all data
disp('MLCC + ReCat is being performed');
imageFFT_mlcc = mlcc_compress(imageFFT, Dp);

% ReCat over MLCC compressed coils
% interpolation kernel estimation
disp('Estimating interpolation kernel - this may take a few minutes');
% crop the calibration region
kCalib = crop(imageFFT_mlcc, [calibSize,N,Dp]);
% estimate the kernel over aggregated tensor
kernel = zeros([kernelSize, N, Dp, N, Dp]);
YtY = data2YtY(kCalib, kernelSize);
for n=1:N
	for d=1:Dp
		kernel(:,:,:,:,n,d) = calibrate(YtY, kernelSize, N, Dp, n, d, beta);
	end
end
% LSQR reconstruction
disp('Performing LSQR reconstruction');
% perform the optimization
res = recat_optimize(imageFFT_mlcc, kernel, nIter, lambda, wavWeight);
% transform to image domain
res_ims = ifft2c( res);
% P-norm combination
result_mlcc = sos(sos(res_ims,4,p_coils),3,p_acq);


%% Normalization
% normalize both reference image and reconstructed images to 0-1.
originalImage = normalize(originalImage);
result_gcc = normalize(result_gcc);
result_mlcc = normalize(result_mlcc);

%% Display
figure; imshow(originalImage);
title('Fully Sampled - Not Compressed Image');
figure; imshow(result_gcc);
title('GCC + ReCat Reconstructed Image');
figure; imshow(result_mlcc);
title('MLCC + ReCat Reconstructed Image');

%% PSNR evaluation
fprintf('GCC + ReCat PSNR: %.2f\n', psnr(result_gcc, originalImage));
fprintf('MLCC + ReCat PSNR: %.2f\n', psnr(result_mlcc, originalImage));
