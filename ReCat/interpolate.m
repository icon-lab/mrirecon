function res = interpolate(x, imkernel, adjoint)
% function to perform interpolation and its adjoint
% NOTE: This function requires SPIRiT library.
%
% res = interpolate(x, imkernel, adjoint)
%
%   Inputs:
%           x        - 4D data array with size [sx, sy, N, D] where
%                      N is the number of phase-cycles, D is the coil
%				       number and the images are sx-by-sy
%           imkernel - interpolation kernel (in image domain)
%			adjoint  - adjoint flag (boolean variable)
%
%   output:
%           res 	- 4D interpolated k-space tensor
%
% (based on SPIRiT operator of SPIRiT library)

	[sx, sy, N, D] = size(x);
	res = zeros(size(x));
	imx = ifft2c(x);
	
	if ~adjoint
        for n=1:N
			for d=1:D
				tmpk = imkernel(:,:,:,:,n,d);
				res(:,:,n,d) = sum(sum(tmpk.*imx,4),3);
			end
        end
	else
        for n=1:N
			for d=1:D
				tmpk = squeeze(conj(imkernel(:,:,n,d,:,:)));
				res(:,:,n,d) = sum(sum(tmpk.*imx,4),3);
			end
        end
	end
	res = fft2c(res)-x;
end