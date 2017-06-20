function res = mlcc_compress( imfft, vcNo )
% Multi-Linear Coil Compression
%
% function to perform coil compression based on MLCC model over a single
% slice acquisition.
% NOTE: This function requires TensorLab package.
%
% res = gcc_compress( imfft, slwin, vcNo )
%
%   Inputs:
%           imfft - 4D array with size [sx, sy, N, D] where N is the number
%                   of phase cycles, D is the coil number and the images
%                   are sx-by-sy
%           vcNo  - number of virtual coils
%
%   output:
%           res - 4D array of compressed coils with size [sx, sy, N, vcNo]
%
%   (c) 2017 Erdem Bıyık

    coilsDim = 4;
    [U, ~, ~] = mlsvd(imfft);
    Ur = U{coilsDim}(:,1:vcNo)';
	res = tmprod(imfft, Ur, coilsDim);
end



