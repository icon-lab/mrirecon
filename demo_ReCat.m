%% Reconstruction by Calibration over Tensors (ReCat) Demo
% This is a demo that shows an example implementation of ReCat algorithm
% over bSSFP data with multiple phase-cycled acquisitions and coils. It is
% based on Biyik et. al, "Reconstruction by Calibration over Tensors for
% Multi-Coil Multi-Acquisition Balanced SSFP Imaging". ReCat is a technique
% that leverages the structural information from both phase-cycles and coils.
% An interpolation kernel is estimated over the aggregated data and used to 
% fill missing k-space points.


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
wavWeight = 0.005; % Wavelet threshold term for l1 regularization (0 to disable)
p_acq = 4; % p value for p-norm combination over phase-cycles
p_coils = 2; % p value for p-norm combination over coils

%% Necessary libraries and folders
% SPIRiT V0.3 is required
if exist('crop.m', 'file')==0
	addpath(genpath('ESPIRiT'));
	if exist('crop.m', 'file')==0
		warning('ReCat utilizes SPIRiT V0.3 library. Please download from http://people.eecs.berkeley.edu/~mlustig/Software.html and copy the "ESPIRiT" directory to the same directory as this demo file. ');
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
load('data/invivo_4coil.mat');
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

%% Interpolation kernel estimation
disp('Estimating interpolation kernel - this may take a few minutes');
% crop the calibration region
kCalib = crop(imageFFT, [calibSize,N,D]);
% estimate the kernel over the tensor
kernel = zeros([kernelSize, N, D, N, D]);
YtY = data2YtY(kCalib, kernelSize);
for n=1:N
	for d=1:D
		kernel(:,:,:,:,n,d) = calibrate(YtY, kernelSize, N, D, n, d, beta);
	end
end

%% LSQR reconstruction
disp('Performing LSQR reconstruction');
% perform the optimization
res = recat_optimize(imageFFT, kernel, nIter, lambda, wavWeight);
% transform to image domain
res_ims = ifft2c( res);

%% P-norm combination
result = sos(sos(res_ims,4,p_coils),3,p_acq);

%% Normalization
% normalize both reference image and reconstructed image to 0-1.
originalImage = normalize(originalImage);
result = normalize(result);

%% Display
figure; imshow(originalImage);
title('Fully Sampled Image');
figure; imshow(result);
title('Reconstructed Image');

%% PSNR evaluation
fprintf('PSNR: %.2f\n', psnr(result, originalImage));
