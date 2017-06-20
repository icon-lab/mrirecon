function [res] = CC(DATA, mtx, dim)
% 
% Geometric Decomposition Coil Compression Based on: Zhang et. al MRM 2013;69(2):571-82.
% The function uses compression matrices computed by calcSCCMtx , calcECCMtx 
% calcGCCMtx to project the data onto them. 
%
%
%  Inputs:
%           DATA - a 4D matrix representing [Kx Ky Kz Coils] data
%                   or a 3D matrix representing [Kx Ky COils]
%           mtx  - Aligned compression matrices from calicECCMtx,
%           calcGCCMtx
%           dim  - dimension to perform compression on. If missing then the
%                   function will assume SCC type compression matrix
%
%  Outputs:
%           res - returned data that is rotated to the compressed space.
%                 The first ncc coils are the compressed & aligned virtual coils.
%
%
% See:
%       calcGCCMtx,  calcECCMtx, alignCCMtx
%
%
% Example Use compressing from 8->4 channels along 2nd dimension:
%
%       load brain_8ch
%       dim = 2;
%       calib = crop(DATA,[24,200,8]);
%       gccmtx =  calcGCCMtx(calib,dim,1);
%       gccmtx = gccmtx(:,1:4,:);
%       gccmtx = alignCCMtx(gccmtx);
%       CCDATA = CC(DATA,gccmtx, dim);
%
%
%
% (c) Michael Lustig 2013


% check if 2D or 3D data
if length(size(DATA))==3
    ndims = 2;
    DATA = permute(DATA,[1,2,4,3]);
else
    ndims = 3;
end

ncc = size(mtx,2);




% Assume SCC if no dimension is used. 
if nargin < 3
   if length(size(mtx)) > 2
       error('compression matrix is probably GCC, but no dimension was entered');
   end
   [Nx,Ny,Nz,Nc] = size(DATA);
   DATA = reshape(DATA,[Nx*Ny*Nz,Nc]);
   res = DATA*mtx;
   res = reshape(res,[Nx,Ny,Nz,ncc]);
   % squeeze back to two dimensions
   if ndims ==2
        res = squeeze(res);
   end
   

% GCC based compressions
else 


    % permute data if compression is done on 2 dimension
    % done for simple reusible code.
    if dim==2
        DATA = permute(DATA,[2,1,3,4]);
    end
    
    if dim==3
        DATA = permute(DATA,[3,1,2,4]);
    end
    
    if dim>3 | dim <1
        disp('Error, compression dimension is 1 2 or 3')
    end
    
    % perform IFFT
    [Nx,Ny,Nz,Nc] = size(DATA);
    im = ifftc(DATA,1);
    res = zeros(Nx,Ny,Nz,ncc);
    
    % rotate data by compression matrices
    for n=1:Nx
        tmpc = reshape(squeeze(im(n,:,:,:)),[Ny*Nz,Nc]);
        res(n,:,:,:) = permute(reshape(tmpc*mtx(:,:,n),[Ny,Nz,ncc]),[4,1,2,3]);
    end
    
    % perform fft
    res = fftc(res,1);
    
    % permute if necessary
    if dim==2
        res = permute(res,[2,1,3,4]);
    end
    
    if dim==3
        res = permute(res,[3,1,2,4]);
    end

end

% squeeze back to two dimensions
if ndims ==2
    res = squeeze(res);
end

