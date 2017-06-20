function x = softThresh(y,t)
    % apply joint sparsity soft-thresholding 
    absy = sqrt(sum(abs(y).^2,3));
    unity = y./(repmat(absy,[1,1,size(y,3)])+eps);

    res = absy-t;
    res = (res + abs(res))/2;
    x = unity.*repmat(res,[1,1,size(y,3)]);
end