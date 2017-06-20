function mtx = alignCCMtx(mtx,ncc)
% mtx = alignCCmtx(mtx, [ ncc)
%
% Align coil compression matrices based on nearest spanning vectors in
% subspaces. This is an implementation based on Zhang et. al MRM 2013;69(2):571-82.
% 
% Inputs:
%           mtx - [nc,nc,slices] input compression matrices
%           ncc - number of  compressed virtual coils to align (default
%           all)
% Outputs:
%           mtx - aligned compression matrices. 
%
%
% See:
%       calcGCCMtx, calcECCMtx, CC
%
% (c) Michael Lustig 2013




[sx,sy,nc] = size(mtx);

if nargin < 2
    ncc = sy;
end

if (sx == sy) & (ncc == sx)
    disp('Warning: number of aligned coils is the same as physical coils, this will do nothing. Either crop mtx or align fewer coils');
end





% align everything based on the middle slice. 
n0 = floor(nc/2);
A00 = mtx(:,1:ncc,n0);

% Align backwards to first slice
A0 = A00;
for n = n0-1:-1:1
    A1 = mtx(:,1:ncc,n);
    C = A1'*A0;
    [U,S,V]= svd(C,'econ');
    P = V*U'; 
    mtx(:,1:ncc,n) = A1*P';
    A0 = mtx(:,1:ncc,n);
end

% Align forward to end slice
A0 = A00;
for n = n0+1:nc
    A1 = mtx(:,1:ncc,n);
    C = A1'*A0;
    [U,S,V]= svd(C,'econ');
    P = V*U';
    mtx(:,1:ncc,n) = A1*P';
    A0 = mtx(:,1:ncc,n);
end

