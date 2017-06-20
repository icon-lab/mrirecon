function res = mtimes(a,b)

if isa(a,'Wavelet') == 0
    error('In  A.*B only A can be Wavelet operator');
end

%res = b*0;

S = size(b);
b = reshape(b,S(1),S(2),prod(S(3:end)));


for n=1:size(b,3)

if a.TI==0
    if a.adjoint
        res(:,:,n) = IWT2_PO(real(b(:,:,n)),a.wavScale,a.qmf) + i*IWT2_PO(imag(b(:,:,n)),a.wavScale,a.qmf);
    else
        res(:,:,n) = FWT2_PO(real(b(:,:,n)),a.wavScale,a.qmf) + i* FWT2_PO(imag(b(:,:,n)),a.wavScale,a.qmf);
    end
else
    if a.adjoint
        res(:,:,n) = IWT2_TI(real(b(:,:,n)),a.wavScale,a.qmf) + i*IWT2_TI(imag(b(:,:,n)),a.wavScale,a.qmf);
    else
        res(:,:,n) = FWT2_TI(real(b(:,:,n)),a.wavScale,a.qmf) + i* FWT2_TI(imag(b(:,:,n)),a.wavScale,a.qmf);
    end
end
end

b = reshape(res,[size(res,1),size(res,2),S(3:end)]);


