function y = normalize( x )
    
    y = x - min(x(:));
    y = y/max(x(:));

end