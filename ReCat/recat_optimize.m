function res = recat_optimize(DATA, kernel, nIter, lambda, wavWeight)
% Optimizer function for ReCat algorithm
% NOTE: This function requires SPIRiT library.
%
% res = recat_optimize(DATA, kernel, nIter, lambda)
%
%   Inputs:
%           DATA   - 4D data array with size [sx, sy, N, D] where
%                    N is the number of phase-cycles, D is the coil
%				     number and the images are sx-by-sy
%           kernel - k-space interpolation kernel
%			nIter  - number of LSQR iterations
% 			lambda - l_2 regularization term weight
%
%   output:
%           res 	- 4D tensor computed by the optimization
%
% (based on several functions and operators of SPIRiT library)

    if nargin < 5
        wavWeight = 0;
    end

	[sx, sy, N, D] = size(DATA);
	ssx = 2^ceil(log2(sx)); 
	ssy = 2^ceil(log2(sy));
	ss = max(ssx, ssy);
    W = Wavelet('Daubechies',4,4);
    
	kernel = kernel(end:-1:1,end:-1:1,:,:,:,:);
	imkernel = zeros(sx, sy, N, D, N, D);
	for n=1:N
		for d=1:D
			imkernel(:,:,:,:,n,d) = ifft2c(zpad(kernel(:,:,:,:,n,d)*sqrt(sx*sy), sx, sy, N, D));
		end
    end
	
	idx_acq = find(abs(DATA)>0);
	idx_nacq = find(abs(DATA)==0);
	empNum = length(idx_nacq(:));
	
	target = interpolate(DATA,imkernel,0);
	target = [-target(:); zeros(empNum,1)];
    xn = zeros(empNum,1);
    for i=1:nIter
        [xn, ~] = lsqr(@interp,target,1e-6,1,speye(empNum,empNum),speye(empNum,empNum),xn,imkernel,sx,sy,N,D,idx_nacq,lambda);
        if wavWeight > 0
            X = DATA;
            X(idx_nacq) = xn;
            X = ifft2c(X);
            X = zpad(X,ss,ss,N,D);
            X = reshape(X,[ss,ss,N*D]);
            X = W*(X);
            X = softThresh(X, wavWeight);
            X = W'*(X);
            X = reshape(X,[ss,ss,N,D]);
            X = crop(X,sx,sy,N,D);
            X = fft2c(X);
            xn = X(idx_nacq);
        end
    end
	res = DATA;
	res(idx_nacq) = xn;
	
	function [res,tflag] = interp(x,imkernel,sx,sy,N,D,idx_nacq,lambda,tflag)
		if ~strcmp(tflag,'transp')
			tmpx = zeros(sx,sy,N,D);
			tmpx(idx_nacq) = x;
			res = interpolate(tmpx,imkernel,0);
			res = [res(:) ; sqrt(lambda)*x(:)];
		else
			tmpy = reshape(x(1:sx*sy*N*D),sx,sy,N,D);
			res = interpolate(tmpy,imkernel,1);
			res = res(idx_nacq) + x(sx*sy*N*D+1:end)*sqrt(lambda);
		end
	end
		
end