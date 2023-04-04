%% Apply gaussian filter to image
function I = imgaussian(I,sigma)

    if sigma==0; return; end; 
    

    radius = ceil(3*sigma);
    [x,y]  = ndgrid(-radius:radius,-radius:radius); % kernel coordinates
    h      = exp(-(x.^2 + y.^2)/(2*sigma^2));
    h      = h / sum(h(:));
    
 
    I = imfilter(I,h);

end