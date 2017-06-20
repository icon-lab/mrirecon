function YtY = data2YtY(data, kernelSize)
% function to compute the calibration matrix from calibration data 
%
% YtY = data2YtY(data, kernelSize)
%
%   Inputs:
%           data       - 4D data array with size [sx, sy, N, D]
%                        where N is the number of phase-cycles, D is
%				         the coil number and the images are sx-by-sy
%           kernelSize - the size of the interpolation kernel
%
%   output:
%           YtY 	   - Y*Y' where Y is the aggregated k-space data
%						 and Y' is its transpose
%
% (based on dat2AtA function of SPIRiT library)


	[sx,sy,N,D] = size(data);

	tmp = im2row(data,kernelSize);
	[winSize,winCount,tN,tD] = size(tmp);
	Y = reshape(tmp, [winSize,winCount*tN*tD]);

	YtY = Y'*Y;


end
