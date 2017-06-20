function res = gcc_compress( imfft, slwin, vcNo )
% Geometric Decomposition Coil Compression
%
% function to perform coil compression based on GCC model over a single
% slice acquisition.
% NOTE: This function requires SPIRiT library.
%
% res = gcc_compress( imfft, slwin, vcNo )
%
%   Inputs:
%           imfft - 3D array with size [sx, sy, D] where D is the coil
%                   number and the images are sx-by-sy
%           slwin - sliding window length
%           vcNo  - number of virtual coils
%
%   output:
%           res - 3D array of compressed coils with size [sx, sy, vcNo]
%
% (c) 2017 Erdem Bıyık (based on GCC demo from SPIRiT library)

    readoutDim = 3;
    imfft = permute(imfft,[1 2 4 3]);
    gccmtx = calcGCCMtx(imfft, readoutDim, slwin);
    
    % No alignment step in compression for single slice 2D undersampling case 
    
    res = CC(imfft,gccmtx(:,1:vcNo), readoutDim);
    res = permute(res,[3 1 4 2]);
end

