function kernel = calibrate(YtY, kernelSize, N, D, n, d, beta)
% function to solve Y*t = y for t with Tikhonov regularization
%
% kernel = calibrate(YtY, kernelSize, N, D, n, d, beta)
%
%   Inputs:
%           YtY        - Y*Y' matrix where Y' is the transpose of Y
%           kernelSize - the size of the kernel to be approximated
%           N          - number of phase-cycles
%           D          - number of coils
%           N          - id of the current phase-cycle
%           N          - id of the current coil
%			beta       - Tikhonov regularization parameter
%
%   output:
%           kernel     - 4D array that corresponds to the
%						 interpolation kernel for n-th phase-cycle
% 						 and d-th coil
%
% (based on the calibrate function of SPIRiT library)

	idxy = prod(kernelSize)*((d-1)*N + (n-1) + 0.5) + 0.5;
	idxY = [1:idxy-1, idxy+1:prod(kernelSize)*N*D];

	Yty = YtY(:,idxy);
	Yty = Yty(idxY);
	YtY = YtY(idxY,:);
	YtY = YtY(:,idxY);

	kernel = zeros([kernelSize,N,D]);
	beta = norm(YtY,'fro')/size(YtY,1)*beta;
	kernel(idxY) = inv(YtY + eye(size(YtY))*beta)*Yty;	
	kernel = reshape(kernel, [kernelSize,N,D]);

end












