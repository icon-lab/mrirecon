function mtx = calcGCCMtx(calibDATA,dim,ws)
% mtx = calcGCC(calibDATA,dim,[ws )
%
% Geometric Decomposition Coil Compression Based on: Zhang et. al MRM 2013;69(2):571-82.
% The function computes and returns compression matrices for each position in space
% along a chosen dimension. To compress, only use the Nth first entries of the 2st
% dimension of mtx. Use with alignCCMtx to align and CC to perform compression
%
%
%  Inputs:
%           calibDATA - a 4D matrix representing [Kx Ky Kz Coils] calibration data
%                   or a 3D matrix representing [Kx Ky COils]
%           dim  - dimension to perform compression on
%           ws   - odd number of window size in space for computing the cc matrices
%                   default is ws=1
%
%  Outputs:
%           mtx - the compression matrix.
%
% See:
%       calcECCMtx, ECC, CC, alignCCMtx
%
% (c) Michael Lustig 2013


if nargin < 3
    ws = 1;
end


% check if k-space is 2D or 3D 
if length(size(calibDATA))==3
    ndims = 2;
    calibDATA = permute(calibDATA,[1,2,4,3]);
else
    ndims = 3;
end


% round sliding window size to nearest odd number
ws = floor(ws/2)*2 + 1;


% permute data if compression is done on 2nd or 3rd dimensions
% done for simple reusible code.
if dim==2
    calibDATA = permute(calibDATA,[2,1,3,4]);
end

if dim==3
    calibDATA = permute(calibDATA,[3,1,2,4]);
end

if dim>3 | dim <1
    disp('Error, compression dimension is 1 2 or 3')
end

% perform IFFT
[Nx,Ny,Nz,Nc] = size(calibDATA);
im = ifftc(calibDATA,1);
res = zeros(Nx,Ny,Nz,Nc);


% check if there's enough data for subspace calculation
if ws*Ny*Nz < 3* Nc
    disp('Warning: ratio between data in each slice and number of channel is less than 3 -- noise could bias results. You should increas ws')
end


% Calculate compression matrices for each position in the readout
% over a sliding window of size ws


mtx = zeros(Nc,min(Nc,ws*Ny*Nz),Nx);
zpim = zpad(im,[Nx + ws-1,Ny,Nz,Nc]);
for n = [1:Nx]
    %tmpc = squeeze(im(n,:,:));
    tmpc = reshape(zpim(n:n+ws-1,:,:,:),ws*Ny*Nz,Nc);
    [U,S,V] = svd(tmpc,'econ');
    mtx(:,:,n) = V;
end

